classdef SVOR
    %SVOR Summary of this class goes here
    %Detailed explanation goes here
    
    properties
        alpha=[];   %dual variables
        b=0;        %offset
        pms=0;       %one parameters C
        mu=0;        % mu
        size_subclass=0;         %the size of each class
        active_size=0;
        iteration=0;
        KKT=0;        % kkt condistion: 1, yes; 0 no
    end
    methods(Static = true)
        %-----------  add J
        function obj=SVOR(X,Y,C)
            global  x y KTYPE KSCALE initUbIn alpha_home j;
            %-----------预处理
            [origin.y,sort_index]=sort(Y);
            unique_label=unique(origin.y);
            unique_num=hist(origin.y,unique_label);
            
            l=2*length(Y)-unique_num(1)-unique_num(end);
            copy.x=[];
            copy.y=[];
            copy.j=[];
            sum_x=0;
            for i=1:length(unique_label)-1
                copy.x=[copy.x;X([find(Y==i);find(Y==i+1)],:)];
%                 copy.x=[copy.x;X(sort_index(sum_x+1:sum_x+unique_num(i)+unique_num(i+1)),:)];
                copy.y=[copy.y;-ones(unique_num(i),1);ones(unique_num(i+1),1)];
                copy.j=[copy.j;ones(unique_num(i)+unique_num(i+1),1)*i];
                sum_x=sum_x+unique_num(i);
            end
            
            alpha0=zeros(l,1);
            clear unique_label  i l origin sort_index sum_x
            %-----------预处理结束
            x=copy.x;
            y=copy.y;
            j=copy.j;
            clear copy;
%             KTYPE = ktype;
%             KSCALE =kscale;
            initUbIn=C;  %set the upper bound
            
            [mu,b,active_size,iteration,kkt]=SVOR.quadsmo(alpha0);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.mu=mu;
            obj.pms=C;
            obj.alpha=alpha_home;
            obj.b=b;
            obj.active_size=active_size';
            obj.iteration=iteration;
            obj.size_subclass=unique_num;
            obj.KKT=kkt;

        end
        function k = Kernel(x, y)%the kernel function
            % function k = kernel(x, y);
            %
            %	x: (Lx,N) with Lx: number of points; N: dimension
            %	y: (Ly,N) with Ly: number of points
            %	k: (Lx,Ly)
            %
            %	KTYPE = 1:      linear kernel:      x*y'
            %	KTYPE = 2,3,4:  polynomial kernel:  (x*y'*KSCALE+1)^KTYPE
            %	KTYPE = 5:      sigmoidal kernel:   tanh(x*y'*KSCALE)
            %	KTYPE = 6:	gaussian kernel with variance 1/(2*KSCALE)
            %
            %       assumes that x and y are in the range [-1:+1]/KSCALE (for KTYPE<6)
            
            global KTYPE
            global KSCALE
            
            k = x*y';
            if KTYPE == 1				% linear
                % take as is
            elseif KTYPE <= 4			% polynomial
                k = (k*KSCALE+1).^KTYPE;
            elseif KTYPE == 5			% sigmoidal
                k = tanh(k*KSCALE);
            elseif KTYPE == 6			% gaussian
                [Lx,~] = size(x);       % the number of x rows
                [Ly,~] = size(y);
                k = 2*k;
                k = k-sum(x.^2,2)*ones(1,Ly);   %sum(A,2) means compute the sum of the elements in each row
                k = k-ones(Lx,1)*sum(y.^2,2)';
                k = exp(k*KSCALE);
            end
        end
        function [mu,b,a_size,index_iteration,outflag] = quadsmo(alpha0)
            %-----------
            %-------------add j mu j_total b
            global fake_zero x y max_iteration_smo epsilon sample_length index_home shrink_state...
                initUbIn flag_up flag_low  cache_size  regrad_count cache_memory alpha_home j mu j_total B_up B_low
            
            max_iteration_smo=1000000;
            fake_zero=10^-10;
            epsilon=10^-10;
            sample_length=length(y);%n
            cache_memory=40;%MB
            
            counter=min(sample_length,1000)+1;
            index_home=(1:sample_length)';
            alpha_home=alpha0;%n*1  initialize the alpha
            flag_up=((alpha0>fake_zero) & y==-1) | ((alpha0<initUbIn -fake_zero) & y==1);%I_up
            flag_low=((alpha0>fake_zero) & y==1) | ((alpha0<initUbIn -fake_zero) & y==-1);%I_low
            %-------------
            %-------------
            index_iteration=1;
            b=0;
            stopnum=0;
            Q_Li=zeros(sample_length,1);
            regrad_count=0;
            
            used_cache=0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            shrink_state=0;
            %-----------
            %初始化判别函数 -y_i*grad(f)_i & gi_bar
            local_gi=y;
            gi_bar=zeros(size(y));
            active_size=sample_length;
            %-----------
            %-----------
            
            %初始化活动集的大小和索引
            shrink_count=0;
            active_index=(1:sample_length)';%假定的
            a_size=active_size;%假定的
            local_j=j(active_index);
            mu=zeros(j(end)-j(1),1);
            j_total=j(end)-j(1)+1;
            B_up=zeros(j(end)-j(1),1);
            B_low=B_up;
            local_flag_up=flag_up;
            local_flag_low=flag_low;
            local_x=x;
            local_y=y;
            local_alpha=alpha_home;
            %-----------
            %%%%%%%%%%%%%%初始化缓存
            cache_size=floor(2^17*cache_memory/active_size);
            Cache_ind=repmat(int32(0),cache_size,1);
            Cache_iter=repmat(int32(0),cache_size,1);
            main_problem_fail=0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            while index_iteration < max_iteration_smo
%                 counter=counter-1;
%                 if(counter==0)      %%%%%缩减
%                     counter=min(sample_length,1000)+1;
%                     if main_problem_fail==1
%                         main_problem_fail=0;
%                         flag_up(active_index)=local_flag_up;
%                         flag_low(active_index)=local_flag_low;
%                         alpha_home(active_index)=local_alpha;
%                         [leave_index,active_index,active_size,local_gi]=SVOR.DoShrinking(gi_home,gi_bar,flag_up,flag_low,alpha_home,y,index_home,j);
%                         shrink_count=shrink_count+1;
%                         a_size(shrink_count)=active_size;
%                         other_active_index=(1:sample_length);
%                         other_active_index(active_index)=[];
%                         
%                         clear Cache;
%                         clear Cache_ind;
%                         clear Cache_iter;
%                         cache_size=floor(2^17*cache_memory/active_size);
%                         Cache_ind=repmat(int32(0),cache_size,1);
%                         Cache_iter=repmat(int32(0),cache_size,1);
%                         used_cache=0;
%                     else
%                         %%%保存之前的结果
%                         alpha_home(active_index)=local_alpha;
%                         flag_up(active_index)=local_flag_up;
%                         flag_low(active_index)=local_flag_low;
%                         [leave_index,active_index,active_size,local_gi]=SVOR.DoShrinking(local_gi,gi_bar,local_flag_up,local_flag_low,local_alpha,local_y,active_index,local_j);
%                         shrink_count=shrink_count+1;
%                         a_size(shrink_count)=active_size;
%                         other_active_index=(1:sample_length);
%                         other_active_index(active_index)=[];
%                         
%                         if length(leave_index)==sample_length
%                             clear Cache;
%                             clear Cache_ind;
%                             clear Cache_iter;
%                             cache_size=floor(2^17*cache_memory/active_size);
%                             Cache_ind=repmat(int32(0),cache_size,1);
%                             Cache_iter=repmat(int32(0),cache_size,1);
%                             used_cache=0;
%                         else
%                             for i=1:used_cache
%                                 Cache{i}=Cache{i}(leave_index);
%                             end
%                         end
%                     end
%                     %%更新当前活动集
%                     local_flag_up=flag_up(active_index);
%                     local_flag_low=flag_low(active_index);
%                     local_x=x(active_index,:);
%                     local_y=y(active_index);
%                     local_alpha=alpha_home(active_index);
%                     %-----------
%                     local_j=j(active_index);
%                     %-----------
%                     %%%%%%更新local_gi
%                     local_gi=local_gi(leave_index);
%                 end    %%%%%%%%缩减完成
                %_______________%%%选择两个变量
%                 if index_iteration==4420
%                     a=1;
%                 end
                [max_delta_B,local_index,select_j]=SVOR.SelectTwoAlpha(local_gi,local_j,local_flag_up,local_flag_low);
%                fprintf('%d\n',index_iteration);
                %---------------
                %更新两个alpha
                %%%查找缓存块
                cache_index1=find(Cache_ind==active_index(local_index(1)),1);
                cache_index2=find(Cache_ind==active_index(local_index(2)),1);
                
                if isempty(cache_index1) && isempty(cache_index2)
                    if used_cache==cache_size
                        initQAB=SVOR.Kernel(local_x,(local_x(local_index,:))).*(local_y*(local_y(local_index))');%n*2
                        [~,replace_index]=sort(Cache_iter);
                        Cache{replace_index(1)}=initQAB(:,1);
                        Cache{replace_index(2)}=initQAB(:,2);
                        Cache_ind(replace_index([1,2]))=active_index(local_index);
                        Cache_iter(replace_index([1,2]))=index_iteration;
                    elseif used_cache==cache_size-1
                        initQAB=SVOR.Kernel(local_x,(local_x(local_index(1),:))).*(local_y*(local_y(local_index(1)))');
                        Cache{used_cache+1}=initQAB(:,1);
                        Cache_ind(used_cache+1)=active_index(local_index(1));
                        Cache_iter(used_cache+1)=index_iteration;
                        used_cache=used_cache+1;
                        initQAB(:,2)=SVOR.Kernel(local_x,(local_x(local_index(2),:))).*(local_y*(local_y(local_index(2)))');
                        [~,replace_index]=min(Cache_iter);
                        Cache{replace_index}=initQAB(:,2);
                        Cache_ind(replace_index)=active_index(local_index(2));
                        Cache_iter(replace_index)=index_iteration;
                    else
                        initQAB=SVOR.Kernel(local_x,(local_x(local_index,:))).*(local_y*(local_y(local_index))');%n*2
                        
                        Cache{used_cache+1}=initQAB(:,1);
                        Cache{used_cache+2}=initQAB(:,2);
                        Cache_ind([used_cache+1,used_cache+2])=active_index(local_index);
                        Cache_iter([used_cache+1,used_cache+2])=index_iteration;
                        used_cache=used_cache+2;
                    end
                else
                    if cache_index1
                        initQAB=Cache{cache_index1};
                        Cache_iter(cache_index1)=index_iteration;
                    else
                        if used_cache==cache_size
                            initQAB=SVOR.Kernel(local_x,(local_x(local_index(1),:))).*(local_y*(local_y(local_index(1)))');
                            [~,replace_index]=min(Cache_iter);
                            if replace_index==cache_index2
                                [~,replace_index]=sort(Cache_iter);
                                Cache{replace_index(2)}=initQAB(:,1);
                                Cache_ind(replace_index(2))=active_index(local_index(1));
                                Cache_iter(replace_index(2))=index_iteration;
                                
                            else
                                Cache{replace_index}=initQAB(:,1);
                                Cache_ind(replace_index)=active_index(local_index(1));
                                Cache_iter(replace_index)=index_iteration;
                            end
                        else
                            initQAB=SVOR.Kernel(local_x,(local_x(local_index(1),:))).*(local_y*(local_y(local_index(1)))');
                            Cache{used_cache+1}=initQAB(:,1);
                            Cache_ind(used_cache+1)=active_index(local_index(1));
                            Cache_iter(used_cache+1)=index_iteration;
                            used_cache=used_cache+1;
                            
                        end
                    end
                    if cache_index2
                        initQAB(:,2)=Cache{cache_index2};
                        Cache_iter(cache_index2)=index_iteration;
                    else
                        if used_cache==cache_size
                            initQAB(:,2)=SVOR.Kernel(local_x,(local_x(local_index(2),:))).*(local_y*(local_y(local_index(2)))');
                            [~,replace_index]=min(Cache_iter);
                            Cache{replace_index}=initQAB(:,2);
                            Cache_ind(replace_index)=active_index(local_index(2));
                            Cache_iter(replace_index)=index_iteration;
                        else
                            initQAB(:,2)=SVOR.Kernel(local_x,(local_x(local_index(2),:))).*(local_y*(local_y(local_index(2)))');
                            Cache{used_cache+1}=initQAB(:,2);
                            Cache_ind(used_cache+1)=active_index(local_index(2));
                            Cache_iter(used_cache+1)=index_iteration;
                            used_cache=used_cache+1;
                        end
                    end
                end
                % %_______________计算QAB
                % %_______________
                initQBB=initQAB(local_index,:);
                
                initf=-local_y(local_index).*local_gi(local_index)-initQBB*local_alpha(local_index);
                
                sum_two_alpha=sum(local_alpha(local_index).*local_y(local_index)); %y1alpha1+y2alpha=zeta(constant)
                old_alpha=local_alpha(local_index);
                
                [initAlpha] = SVOR.OneSmo(initQBB,initf,sum_two_alpha,local_index,initUbIn,local_y);
                
%                 fval0=0.5*local_alpha(local_index)'*initQBB*local_alpha(local_index)+initf'*local_alpha(local_index);%the old value of two original alpha               
%                 [initAlpha] = SVOR.OneSmo(initQBB,initf,sum_two_alpha,local_index,initUbIn,local_y);
%                 %%%%能量函数 
%                 fval1=0.5*initAlpha'*initQBB*initAlpha+initf'*initAlpha;  %the new value of the two updated alpha
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 if fval0-fval1<-fake_zero   %error
%                     error('myquadprog do not convergence!! \n')
%                 end
%                 fprintf('%e\n ',fval0-fval1);           
                
                
                new_alpha=initAlpha;
                
                local_alpha(local_index)=initAlpha; %update the two new alpha
                
                
                %%%%%update the new local_flag_up and local_flag_low
                local_flag_up(local_index)=(((initAlpha>fake_zero) & local_y(local_index)==-1) | ((initAlpha<initUbIn -fake_zero) & local_y(local_index)==1));
                local_flag_low(local_index)=((initAlpha>fake_zero) & local_y(local_index)==1) | ((initAlpha<initUbIn -fake_zero) & local_y(local_index)==-1);
                
                
                local_gi=-local_gi.*local_y;
                local_gi=local_gi+initQAB*(new_alpha-old_alpha);
                local_gi=-local_gi.*local_y;
                %_______________update mu
                delta_alpha=(new_alpha(1)-old_alpha(1)).*local_y(local_index(1));
                if ~(select_j(1)==select_j(2))  %注意：这里的mu从下标2开始到j_total
                    if select_j(1)<select_j(2)
                        mu((min(select_j)):(max(select_j)-1))= mu((min(select_j)):(max(select_j)-1))-delta_alpha;
                    else
                        mu((min(select_j)):(max(select_j)-1))=mu((min(select_j)):(max(select_j)-1))+delta_alpha;
                    end
                end
                %_______________
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%更新gi_bar
                
                uij=(abs(old_alpha-initUbIn)<fake_zero);
                update_uij=(abs(new_alpha-initUbIn)<fake_zero);
                
                if active_size==sample_length
                    if ~(uij(1)==update_uij(1))
                        Q_Li=initQAB(:,1);
                        if uij(1)
                            gi_bar=gi_bar-old_alpha(1)*Q_Li;
                        else
                            gi_bar=gi_bar+new_alpha(1)*Q_Li;
                        end
                    end
                    if ~(uij(2)==update_uij(2))
                        Q_Li=initQAB(:,2);
                        if uij(2)
                            gi_bar=gi_bar-old_alpha(2)*Q_Li;
                        else
                            gi_bar=gi_bar+new_alpha(2)*Q_Li;
                        end
                    end
                else
                    if ~(uij(1)==update_uij(1))
                        Q_Li(active_index,:)=initQAB(:,1);
                        Q_Li(other_active_index,:)=SVOR.Kernel(x(other_active_index,:),(x(active_index(local_index(1)),:))).*(y(other_active_index)*(y(active_index(local_index(1))))');
                        if uij(1)
                            gi_bar=gi_bar-old_alpha(1)*Q_Li;
                        else
                            gi_bar=gi_bar+new_alpha(1)*Q_Li;
                        end
                    end
                    if ~(uij(2)==update_uij(2))
                        Q_Li(active_index,:)=initQAB(:,2);
                        Q_Li(other_active_index,:)=SVOR.Kernel(x(other_active_index,:),(x(active_index(local_index(2)),:))).*(y(other_active_index)*(y(active_index(local_index(2))))');
                        if uij(2)
                            gi_bar=gi_bar-old_alpha(2)*Q_Li;
                        else
                            gi_bar=gi_bar+new_alpha(2)*Q_Li;
                        end
                    end
                end
                
                %%%%%%判断
                %----------------
%                 fprintf('sub_mM:%e\n ',max_delta_B);
                if max_delta_B<=epsilon   %stopping condition (subproblem)
                    stopnum=stopnum+1;
                    fprintf('sub_mM:%e\n ',max_delta_B);
                    alpha_home(active_index)=local_alpha;
                    flag_up(active_index)=local_flag_up;
                    flag_low(active_index)=local_flag_low;
                    if active_size==sample_length;
                        gi_home=local_gi;
                    else
                        [gi_home]=SVOR.Reconstruct_gradient(local_gi,gi_bar,active_index);
                    end
                    [max_delta_B,~,~]=SVOR.SelectTwoAlpha(gi_home,j,flag_up,flag_low);
                    if (max_delta_B<=epsilon)
                        b=(B_up+B_low)/2;
                        fprintf('main_mM:%e\n ',max_delta_B);
                        break;
                    else
                        main_problem_fail=1;
                        counter=1;
                    end
                end
                
                %----------------
                %%%%%%%%%%%%%%%%%%%%%%%%
                index_iteration = index_iteration+1;
            end
            
            fprintf('%d\n',index_iteration);
            if index_iteration >= max_iteration_smo
                fprintf('index_iteration >= max_iteration_smo\n');
            end
%             %%%%%%%%%%%%%%%%%%%%%%%%%check the KKT conditions
            local_obj.alpha=alpha_home;
            local_obj.b=b;
            local_obj.mu=mu;
            local_type=unique(j);
            every_j_num=hist(j,local_type);
            outflag=SVOR.TestKKT(local_obj,initUbIn,every_j_num);
            if outflag==1
                fprintf('meet the KKT conditions\n');
            else
                fprintf('don not meet the KKT conditions\n');
            end
        end
        %---------------
        function [max_delta_B,local_index,select_j]=SelectTwoAlpha(local_gi,local_j,local_flag_up,local_flag_low)
            global mu j_total B_up B_low
            local_j_type=unique(local_j);
            every_j_num=hist(local_j,local_j_type);
            sum=0;
            
            b_up=zeros(j_total,1);
            b_low=b_up;
            
            B_bar_low=zeros(j_total,1);
            B_bar_up=zeros(j_total,1);
            for j=1:j_total
                local_gi_j=local_gi(sum+1:sum+every_j_num(j));
                local_up_j=local_flag_up(sum+1:sum+every_j_num(j));
                local_low_j=local_flag_low(sum+1:sum+every_j_num(j));
                
                index_j=sum+1:sum+every_j_num(j);
                index_up_j=index_j(local_up_j);
                index_low_j=index_j(local_low_j);
                if any(local_up_j)
                [b_up(j),ind_b_up_j]=max(local_gi_j(local_up_j));
                index_b_up_j(j)=index_up_j(ind_b_up_j);
                else
                b_up(j)=-inf;
%                   index_b_up_j(j)=  1;
                end
                if any(local_low_j)
                [b_low(j),ind_b_low_j]=min(local_gi_j(local_low_j));
                index_b_low_j(j)=index_low_j(ind_b_low_j);
                else
                    b_low(j)=inf;
                end
                sum=sum+every_j_num(j);
            end
            for j=1:j_total
                [B_bar_low(j),ind_B_bar_low]=min(b_low(1:j));
                [B_bar_up(j),ind_B_bar_up]=max(b_up(j:end));
                index_B_bar_low(j)=index_b_low_j(ind_B_bar_low);
                index_B_bar_up(j)=index_b_up_j(ind_B_bar_up+j-1);
            end
            B_low=B_bar_low;
            B_up=B_bar_up;
            alpha_low=index_B_bar_low;
            alpha_up=index_B_bar_up;
            for j=1:j_total-1
                if ~(j==j_total-2)
                    if mu(j)>0        %这里的mu从下标2开始到j_total，不用+1
                        B_low(j)=B_bar_low(j+1);
                        alpha_low(j)=index_B_bar_low(j+1);
                    end
                end
                if ~(j==1)
                    if mu(j-1)>0
                        B_up(j)=B_bar_up(j-1);
                        alpha_up(j)=index_B_bar_up(j-1);
                    end
                end
            end
            [max_delta_B,j_index]=max(B_up-B_low);
            local_index=[alpha_up(j_index),alpha_low(j_index)];
            select_j=[local_j(local_index(1)),local_j(local_index(2))];
            
            
        end
        %---------------
        function [gi]=Reconstruct_gradient(local_gi,gi_bar,active_index)
            global x y initUbIn fake_zero sample_length alpha_home regrad_count;
            
            other_active_index=(1:sample_length);
            other_active_index(active_index)=[];
            gi=zeros(sample_length,1);
            gi(active_index,:)=-y(active_index).*local_gi;
            gi(other_active_index,:)=gi_bar(other_active_index)-1;  %初始化非活动集的gi
            free_flag=abs(alpha_home(active_index)-initUbIn)>fake_zero & alpha_home(active_index)>fake_zero;
            free_index=active_index(free_flag);
            
            reSize=100000;
            if length(other_active_index)>reSize
                reGroup=length(other_active_index)/reSize;
                for i=1:reGroup
                    QPiAF=SVOR.Kernel(x(other_active_index(((i-1)*reSize+1):i*reSize),:),x(free_index,:)).*(y(other_active_index(((i-1)*reSize+1):i*reSize))*(y(free_index))');
                    gi(other_active_index(((i-1)*reSize+1):i*reSize))=gi(other_active_index(((i-1)*reSize+1):i*reSize))+QPiAF*alpha_home(free_index);
                end
                if i*reSize<length(other_active_index)
                    QPiAF=SVOR.Kernel(x(other_active_index(i*reSize+1:end),:),x(free_index,:)).*(y(other_active_index(i*reSize+1:end))*(y(free_index))');
                    gi(other_active_index(i*reSize+1:end))=gi(other_active_index(i*reSize+1:end))+QPiAF*alpha_home(free_index);
                end
            else
                QNAf=SVOR.Kernel(x(other_active_index,:),(x(free_index,:))).*(y(other_active_index)*(y(free_index))');
                gi(other_active_index)=gi(other_active_index)+QNAf*alpha_home(active_index(free_flag));
            end
            gi=-y.*gi;
            regrad_count=regrad_count+1;
        end
        function [lea_index,act_index,active_size,local_gi]=DoShrinking(local_gi,gi_bar,lcoal_flag_up,local_flag_low,local_alpha,local_y,active_index,local_j)
            global y fake_zero shrink_state epsilon initUbIn alpha_home flag_up flag_low index_home j_total B_up B_low
            %---------------更改的shrink
            local_j_type=unique(local_j);
            every_j_num=hist(local_j,local_j_type);
            sum=0;
            act_index=[];
            lea_index=[];
            %             if(shrink_state==0) && (max_value-min_value)<=10*epsilon
            %                 shrink_state=1;
            %                 %%%%
            %                 [gi_home]=SVOR.Reconstruct_gradient(local_gi,gi_bar,active_index);
            %                 local_gi=gi_home;
            %                 %%% 最大违反对 这里需要重新求吗？
            %                 [max_value]=max(gi_home(flag_up));
            %                 [min_value]=min(gi_home(flag_low));
            %                 flag_low_bound=(alpha_home>initUbIn-fake_zero & y==1) | (alpha_home<fake_zero & y==-1);
            %                 flag_up_bound=(alpha_home<fake_zero & y==1) | (alpha_home>initUbIn-fake_zero & y==-1);
            %                 leave_index=(gi_home>max_value+fake_zero & flag_low_bound)...
            %                     |(gi_home<min_value-fake_zero & flag_up_bound);
            %                 leave_index=~leave_index;
            %                 active_index=index_home(leave_index);
            %                 active_size=length(active_index);
            %             else
            flag_low_bound=(local_alpha>initUbIn-fake_zero & local_y==1) | (local_alpha<fake_zero & local_y==-1);
            flag_up_bound=(local_alpha<fake_zero & local_y==1) | (local_alpha>initUbIn-fake_zero & local_y==-1);
            for j=1:j_total
                
                leave_index=(local_gi(sum+1:sum+every_j_num(j))>B_up(j)+fake_zero & flag_low_bound(sum+1:sum+every_j_num(j)))...
                    |(local_gi(sum+1:sum+every_j_num(j))<B_low(j)-fake_zero & flag_up_bound(sum+1:sum+every_j_num(j)));
                leave_index=~leave_index;
                lea_index=[lea_index;leave_index];
                a_index_j=active_index((sum+1:sum+every_j_num(j)));
                act_index=[act_index;a_index_j(leave_index)];
                active_size=length(act_index);
                sum=sum+every_j_num(j);
            end
            %             end
            %---------------
            lea_index=logical(lea_index);
        end
        function [initAlpha] = OneSmo(initQ,initf,sum_two_alpha,two_index,ub,local_y) %update the two alpha
            first_index=two_index(1);
            second_index=two_index(2);
            y1=local_y(first_index);
            y2=local_y(second_index);
            yy=y1*y2;
            %convert to quadratic equation of one variable
            a=(initQ(1,1)+initQ(2,2)-2*yy*initQ(1,2))/2;
            b=-sum_two_alpha*y1*initQ(2,2)+initQ(1,2)*y2*sum_two_alpha+initf'*[1;-yy];
            % c=0.5*initQ(2,2)*sum_two_alpha*sum_two_alpha + initf(2)*sum_two_alpha*y2;
            if yy==1    %y1=y2  y1*alpha1+y2*alpha2=zeta --y1(y1*alpha1+y2*alpha)=alpha1+alpha2
                left_bound=max(0,y1*sum_two_alpha-ub);     %L
                right_bound=min(ub,y1*sum_two_alpha);      %H
            else    %y1<>y2  y1(y1*alpha1+y2*alpha)=alpha1-alpha2
                left_bound=max(0,y1*sum_two_alpha);
                right_bound=min(ub,ub+y1*sum_two_alpha);
            end
            opitimal_alpha1=-b/(2*a);
            if opitimal_alpha1 > right_bound
                opitimal_alpha1=right_bound;
            end
            if opitimal_alpha1 < left_bound
                opitimal_alpha1=left_bound;
            end
            initAlpha=[opitimal_alpha1; y2*sum_two_alpha-yy*opitimal_alpha1];
        end
        function out_KKT=TestKKT(os,C,every_j_num)
            %---------------add i_total
            global fake_zero2 x y  j_total
            %---------------
            fake_zero2=10^-3;
            %             out_KKT=0;
            %             local_g=SVOR.CaculateF(os);
            alpha=os.alpha;
            %---------------check KKT : three steps
            out_KKT=1;
            check_mu=os.mu;
            check_b=os.b;
            check_mu=[0;check_mu;0];
            sum=0;
            %%%核函数
            Q=SVOR.Kernel(x,x).*(y*y');
            for j=1:j_total
                %%%check every group alpha'*y_i
                y_alpha_sum=alpha(sum+1:sum+every_j_num(j))'*y(sum+1:sum+every_j_num(j));
                if abs(y_alpha_sum-check_mu(j)+check_mu(j+1))>fake_zero2
                    out_KKT=0;
                    break;
                end
                %%%check every mu
                if ~(j==1)
                    if check_mu(j)==0 && check_b(j-1)-check_b(j)<fake_zero2
                        out_KKT=0;
                        break;
                    end
                    if check_mu(j)>0 && abs(check_b(j-1)-check_b(j))>fake_zero2
                        out_KKT=0;
                        break;
                    end
                end
                
                %%%%calculate local_g_j
                local_Q=[y(sum+1:sum+every_j_num(j)) Q(sum+1:sum+every_j_num(j),:)];
                local_alpha=[os.b(j);os.alpha];
                local_f=local_Q*local_alpha;
                local_g=local_f-1;
                %%%check every alpha
                local_tmp=SVOR.SubKKT(alpha(sum+1:sum+every_j_num(j)),local_g,C); %check every alpha
                if ~(local_tmp==1)
                    out_KKT=0;
                    break;
                end
                sum=sum+every_j_num(j);
            end
            
            %---------------
        end
        function out_KKT=SubKKT(alpha,local_g,C)
            global fake_zero2
            out_KKT=1;
            sum_length=length(alpha);
            for i=1:sum_length
                if alpha(i)<fake_zero2
                    if local_g(i)<-fake_zero2
                        out_KKT=0;
                        break;
                    end
                else
                    if alpha(i)>C-fake_zero2
                        if local_g(i)>fake_zero2
                            out_KKT=0;
                            break;
                        end
                    else
                        if local_g(i)>fake_zero2 || local_g(i)<-fake_zero2
                            out_KKT=0;
                            break;
                        end
                    end
                end
            end
        end
    end
end

