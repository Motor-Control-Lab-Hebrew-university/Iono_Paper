% classify_units
addpath F:\monkeydata\Penny\analysis\IONO2SU\Classifier\prepare_data
indir= 'F:\monkeydata\Penny\analysis\IONO2SU\create_database4allunits\units2class\';
files= dir([indir]);
c_all=[0 0];
ex_ax = [1 11];
for i=1:length(files)
    if exist([indir files(i).name])==2
        load([indir files(i).name]);
        if ~isempty(unit.Waveform.ax)
            [p2t, waveamp, slope,decay,width] = get_waveform_par(unit.Waveform.ax,unit.Waveform.Trc);
            if ~isempty(find(unit.SCP))
                [~,index]=max([unit.pst2stim.Amp]);
                PSTH=[];
                for ii=1:length(index)
                    if ~isempty(unit.pst2stim(index).PSTH) & ~isempty(find(ismember([unit.pst2stim.resp],[1 3])))
                        PSTH=[PSTH; unit.pst2stim(index).PSTH];
                    end
                end
                if ~isempty(PSTH)
                    ax=unit.pst2stim(index).ax;
                    bix = find(ax<-10);
                    baseline=mean(mean(PSTH(:,bix),2));
                    norm_PSTH= PSTH-baseline;
                    ex_ix= find(ax>ex_ax(1) & ax<ex_ax(2));
                    for exi = 1:size(norm_PSTH,1)
                        index= find(norm_PSTH(exi, ex_ix)>=0);
                        control_ex(exi)=sum(norm_PSTH(exi,ex_ix(index)));
                    end
                    excitation = mean(control_ex);
                    pars.width = width;
                    pars.excitation =excitation;
                    edfiles=find(~unit.SCP & ~unit.HFS);
                    if unit.edname(1)=='m' & unit.sess<30
                        continue
                    else
                        baseline= find_baseline(unit.edname,unit.uid,edfiles);
                        
                        pars.Baseline=baseline;
                        
                        
                        if unit.IS>0.1
                            disp(['classifying ' files(i).name '-'])
                            [c,p]=find_class4unit(pars);
                            unit.clean_class.label= c;
                            unit.clean_class.probability= p;
                            if c==1
                                unit.clean_class.name= 'Pyramidal';
                            else
                                unit.clean_class.name= 'PV';
                            end
                            disp(['assigned ' unit.clean_class.name '; score=' num2str(p) '.. and saving..'])
                        else
                            unit.clean_class.label= nan;
                            unit.clean_class.probability= nan;
                        end
                        
                    end
                else
                    unit.clean_class.label= nan;
                    unit.clean_class.probability= nan;
                end
                c_all(c)=c_all(c)+1;
                save([indir files(i).name],'unit');
            end
        end
    end
end