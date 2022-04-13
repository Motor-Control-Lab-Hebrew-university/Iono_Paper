function baseline= find_baseline(edname,uid,files)
ee_index= strfind(edname,'ee');
sessid = str2num(edname(2:ee_index-3));
subid = str2num(edname(ee_index-2:ee_index-1));
date= id2date(sessid, edname(1));
if edname(1) == 'p'
    indir = 'G:\Penny\Ctx_rec\';
elseif edname(1) == 'm'
    indir= 'G:\Menta\CTX_rec\';
else
    indir= 'G:\Carmen\Ctx_rec\';
end
nSpikes=0;
Times=0;
for i=1:length(files)
    if isdir([indir date '\MergedEdFilesSorted\'])
        edfile= load([indir date '\MergedEdFilesSorted\' edname '.' num2str(files(i)) '.mat']);
    else
        edfile= load([indir date '\MergedEdFiles\' edname '.' num2str(files(i)) '.mat']);
    end
    
    if isfield(edfile,['Tspike' num2str(uid)])
       nSpikes=nSpikes + length(edfile.(['Tspike' num2str(uid)]));
       if ~isempty(edfile.TimeEnd)
       Times= Times+ (edfile.TimeEnd-edfile.TimeBegin);
       else %use mat files
           el = floor(uid/100);
           infofile=load([indir date '\Info\' date '_param.mat']);
           filenum=infofile.SESSparam.SubSess(subid).Files(1)+i-1;
           matfile= load([indir date '\MAT\' date sprintf('%03d',filenum) '_wvf.mat'],['CUnit' num2str(el) '_TimeBegin'],['CUnit' num2str(el) '_TimeEnd']);
           Times=Times+ (matfile.(['CUnit' num2str(el) '_TimeEnd'])-matfile.(['CUnit' num2str(el) '_TimeBegin']))*1000;
       end
    end
end
Times=Times/1000;
baseline=nSpikes/Times;
end