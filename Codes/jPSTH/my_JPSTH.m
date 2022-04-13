% res - raw jpsth
%shift_predict - shift predictor using 1 trial
%psth_pred  - psth predictor
% results are in probability of finding a spike
% at the moment, the shift predictor is cyclic
function [ res shift_predict psth_pred surprise_mat std_mat NUM_TRIAL pval_mat] =  ...
    my_JPSTH(cut1, cut2, BIN_SIZE,gauss_filt_std)
% gauss_filt_std is the standard deviation of the gaussian filter for the
% 2d interpulation. By default, or if gauss_filt_std is 0,  no filtering is
% done. 
% pval_mat - matrix with cumulative p-values of having a value in the jpsth
% matrix, with binomial assumption, given the psth-predicted matrix. values
% close to 0 indicate high probability of negative correlation, whereas
% values closer to 1 indicate higher probability of positive correlation
% 
if nargin<3
	gauss_filt_std=0;
end

do_shift=true;
%ignores trials with NaN values, which are considered unsuitable for analysis
trials2take = all(~isnan(cut1))&all(~isnan(cut2)); 
cut1=cut1(:,trials2take);cut2=cut2(:,trials2take); 

NUM_TRIAL = size(cut1,2); %number of trials
CUT_LENGTH = size(cut1,1); %number of samples
BIN_NUM = ceil(CUT_LENGTH/ BIN_SIZE); %number of bins. the last one might not be "broken" (not a whole bin), so it will be thrown at the end.

if mod(CUT_LENGTH,BIN_SIZE)
    bins2throw=1;
else
    bins2throw=0;
end
if NUM_TRIAL <4 %insufficient number of trials make the analysis meaningless.
    
    res = nan(BIN_NUM-bins2throw);
    shift_predict = nan(BIN_NUM-bins2throw);
    psth_pred  = nan(BIN_NUM-bins2throw);
    surprise_mat = nan(BIN_NUM-bins2throw);
    std_mat = nan(BIN_NUM-bins2throw);
    pval_mat = nan(BIN_NUM-bins2throw);
    return
end

% count spikes for JPSTH:

res = zeros(BIN_NUM,BIN_NUM);
shift_predict = zeros(BIN_NUM,BIN_NUM);

NEXT_REF = find(cut2(:,1)); %the first response of cut2 is the first reference
tic
for i=1:NUM_TRIAL %running over the trials
    
    INX1 = find(cut1(:,i)); %running over all spikes in the i'th response of cut1 vector.
    INX2 = NEXT_REF;
    % get next trails
    use_next = i+1;
    if(i == NUM_TRIAL)
        use_next=1;
    end
    NEXT_REF = find(cut2(:,use_next));
    for j=1:length(INX1) %running over all spikes in the i'th response of cut1 vector.
        % update count matrix
        curr_ref = ceil(INX1(j)/BIN_SIZE); %reference index in the binned JPSTH matrix of the j'th spike
        trig = ceil(INX2/BIN_SIZE); %triggered indices in the binned JPSTH matrix

        for m=1:length(trig)
            res(curr_ref, trig(m)) = res(curr_ref, trig(m))+1;
        end
        
        % update shift predictor matrix
        if do_shift
            sp_trig = ceil(NEXT_REF/BIN_SIZE);
            for m=1:length(sp_trig)
                shift_predict(curr_ref, sp_trig(m)) =     shift_predict(curr_ref, sp_trig(m))+1;
            end
        end
                
    end       
end
% t_for=toc

% remove last bin
shift_predict = shift_predict(1:end-bins2throw,1:end-bins2throw);
res = res(1:end-bins2throw,1:end-bins2throw);
    

% get psth predictor
% psth1 = full(sum(cut1',1));
% psth2 = full(sum(cut2',1));
% %bin_psth1 = bin_raster_to_counts(psth1, BIN_SIZE, 0);
%  bin_psth2 = bin_raster_to_counts(psth2, BIN_SIZE, 0);

bin_cut1 = bin_raster_to_counts(cut1', BIN_SIZE, 0);
bin_cut2 = bin_raster_to_counts(cut2', BIN_SIZE, 0);
bin_psth1 = mean(bin_cut1);
bin_psth2 = mean(bin_cut2);
tic
% ress=bin_cut1'*bin_cut2;
% t_mat_mult=toc
std1  = std(bin_cut1);
std2 = std(bin_cut2);
N=2;
std_mat = std1'*std2/(BIN_SIZE/1000)^N;

% psth_pred = bin_psth1'*bin_psth2/(BIN_SIZE/1000)^N; 
psth_pred = bin_psth1'*bin_psth2;

norm_factor = NUM_TRIAL*(BIN_SIZE/1000)^N;
% normalize to probability
% psth_pred = psth_pred/(norm_factor^2);
if(length(psth_pred) ~= length(res))
    psth_pred = psth_pred(1:end-1,1:end-1);
	std_mat = std_mat(1:end-1,1:end-1);
end

if nargout>6
    %     pval_mat = binocdf(res,NUM_TRIAL*BIN_SIZE^2,psth_pred/BIN_SIZE^2);
    pval_mat = poisscdf(res,NUM_TRIAL*psth_pred);
end


LAMDA = psth_pred*norm_factor;
surprise_mat = 1 - poisscdf( res,LAMDA );
surprise_mat  = -log(surprise_mat);
% surprise_mat =[];
shift_predict  = shift_predict / norm_factor;
res = res/norm_factor;
psth_pred = psth_pred/(BIN_SIZE/1000)^N; 
if gauss_filt_std>0
    res  = gauss_filt_2d(res  ,gauss_filt_std);
    shift_predict = gauss_filt_2d(shift_predict ,gauss_filt_std);
    psth_pred = gauss_filt_2d(psth_pred ,gauss_filt_std);
    surprise_mat = gauss_filt_2d(surprise_mat ,gauss_filt_std);
    std_mat = gauss_filt_2d(std_mat ,gauss_filt_std);
end
return;

