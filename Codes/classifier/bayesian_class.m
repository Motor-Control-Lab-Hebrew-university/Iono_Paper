% if ~exist('UNIT_DATA')
%     load ..\..\Clean_data.mat
% end
% clearvars -except UNIT_DATA;
load ..\prepare_data\Classification_parameters_goodunits.mat;
addpath F:\monkeydata\Penny\analysis\ClassificationMethods\RCA\
f=0.6;
iterations=100;
%% prepare parameters
Par2use={'Width' 'Wave_Amp' 'wave_decay' 'AMP' 'Gain' 'Excitation' 'Inhibition' 'Baseline' 'wave_slope' 'Reliability' 'FFresp' 'CVresp'};
Par2use={'Width' 'Baseline' 'Excitation'};
% Par2use={'Width' 'Excitation' 'Baseline'};
Type= [Class_par.Type];
ix1=find(Type==1);
ix2=find(Type==2);
for i=1:length(Par2use)
    ClassPars(:,i)=eval(['[Class_par.' Par2use{i} ']']);
    [~,Weights(i)]=ttest2(ClassPars(ix1,i),ClassPars(ix2,i));
end
ix=find(Vsuccess>0);
ClassPars=ClassPars(ix,:);
Type=Type(ix);
[B,rca,rca_Pars]= RCA(ClassPars,Type);
rca_Pars=ClassPars;
Weights=-log(Weights);
ix = find(~isnan(sum([ClassPars Type(:)],2)));
rca_Pars=rca_Pars(ix,:); Type=Type(ix);
%% classess
mat=zeros(2,3);
ix1=find(Type==1);
ix2=find(Type==2);
for i=1:length(Type)
    if Type(i)==1
        c='Pyr';
    else
        c='FS';
    end
    Class{i}=c;
end
%% class
npars= 3;%[3:size(rca_Pars,2)];%size(rca_Pars,2)];
for nn=npars
    figure
    par2take= combntns([1:size(ClassPars,2)],nn);
    all_ix= [1:size(rca_Pars,1)];
    index= find(ismember(par2take,1));
    par2take=par2take(index,:);
    r=5;
    c=ceil(size(par2take,1)./r);
    for nn1=1:size(par2take,1)
        for i=1:size(rca_Pars,1)
            this_X= rca_Pars(i,par2take(nn1,:));
            for iter=1:iterations
                [~,rix]  = sort(rand(1,length(all_ix)));
                this_ix= all_ix(rix(1:round(f*length(all_ix))));
                this_ix= this_ix(find(~ismember(this_ix,i)));
                this_parameter = rca_Pars(this_ix,par2take(nn1,:));
                this_class = Type(this_ix);
                class_i= Type(i);
                for ii=1:size(this_X,2)
                    %                     [~,p_like]= ttest(this_parameter(find(this_class==1),ii),this_X(ii));
                    %                     p_class= length(find(this_class==1))/length(this_class);
                    %                     [~,p_pre]= ttest(this_parameter(:,ii),this_X(ii));
                    %                     p_pyr_i= (p_like * p_class)/p_pre;
                    %
                    %                     [~,p_like]= ttest(this_parameter(find(this_class==2),ii),this_X(ii));
                    %                     p_class= length(find(this_class==2))/length(this_class);
                    %                     [~,p_pre]= ttest(this_parameter(:,ii),this_X(ii));
                    %                     p_pv_i= (p_like * p_class)/p_pre;
                    p_like= signrank(this_parameter(find(this_class==1),ii)-this_X(ii));
                    p_class= length(find(this_class==1))/length(this_class);
                    p_pre= signrank(this_parameter(:,ii)-this_X(ii));
                    p_pyr_i= (p_like * p_class)/p_pre;
                    
                    p_like= signrank(this_parameter(find(this_class==2),ii)-this_X(ii));
                    p_class= length(find(this_class==2))/length(this_class);
                    p_pre= signrank(this_parameter(:,ii)-this_X(ii));
                    p_pv_i= (p_like * p_class)/p_pre;
                    
                    P_Pyr(iter,ii)=p_pyr_i*Weights(ii);
                    P_PV(iter,ii)=p_pv_i*Weights(ii);
                    
                end
                Mdl = fitcsvm(this_parameter,this_class,'BoxConstraint',10,'KernelFunction','linear');
                svm_class(iter)= (predict(Mdl,this_X));
                %                 Loss(iter)=resubLoss(Mdl);
            end
            p_pv_svm= length(find(svm_class==2))./iterations;
            p_pyr_svm=1-p_pv_svm;
            p_pv_bayes= sum((sum(P_PV,2))>(sum(P_Pyr,2)))./iterations;
            p_pyr_bayes=1-p_pv_bayes;
            classes=[2 1 2 1];
            [~,ix]= max([p_pv_svm p_pyr_svm p_pv_bayes p_pyr_bayes]);
            final_class=classes(ix);
            if p_pv_svm==p_pyr_bayes & p_pyr_svm==p_pv_bayes
                final_class=1;%3;
            end
            
            if final_class==class_i
                mat(class_i,1)=mat(class_i,1)+1;
            elseif final_class~=3
                mat(class_i,2)=mat(class_i,2)+1;
            else
                mat(class_i,3)=mat(class_i,3)+1;
                disp(['Type=' num2str(class_i) ': Width=' num2str(ClassPars(i,1)) '; Baseline=' num2str(ClassPars(i,2)) '; Ex=' num2str(ClassPars(i,3))])
            end
            P_PV=[]; P_Pyr=[];
        end
        subplot(r,c,nn1);
        bar([1 2],mat./[sum(mat,2) sum(mat,2) sum(mat,2)])
        title(Par2use(par2take(nn1,:)));
        mat=zeros(2,2);
    end
%     savefig([num2str(nn) 'pars3']);
%     close
end