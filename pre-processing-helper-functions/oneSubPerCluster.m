%% weighted average across components in each cluster
clusters_to_plot = 3:length(STUDY.cluster);
cluster_lowestIC = STUDY.cluster;
for CL = clusters_to_plot
        rm_comp_ind =[];
        unique_clus_subs = unique(STUDY.cluster(CL).sets);
            for uc = 1:length(unique_clus_subs)
                x = find(STUDY.cluster(CL).sets == unique_clus_subs(uc));
                index = uc;
                if ~isempty(x)
                    if size(x,2)>1 %if subject appears more than once in cluster
                        rm_comp_ind = [rm_comp_ind, x(2:end)];
                    end
                end
            end

new_outlier_comps = STUDY.cluster(CL).comps(rm_comp_ind);
new_outlier_sets =  STUDY.cluster(CL).sets(rm_comp_ind);
cluster_lowestIC(2).comps = [cluster_lowestIC(2).comps, new_outlier_comps]; %CL 2 is outlier comp
cluster_lowestIC(2).sets = [cluster_lowestIC(2).sets, new_outlier_sets];
cluster_lowestIC(CL).comps(rm_comp_ind) = [];
cluster_lowestIC(CL).sets(rm_comp_ind) = [];
end
 disp('Subject IC ERSPs consolidated to')
                        disp('one subject/cluster using lowest IC number')
  STUDY.cluster_og = STUDY.cluster;
  STUDY.cluster_lowestIC = cluster_lowestIC;
  STUDY.cluster =  STUDY.cluster_lowestIC ;
                        
                        
cluster_dist2centroid = STUDY.cluster;
for CL = clusters_to_plot
        rm_comp_ind =[];
        unique_clus_subs = unique(STUDY.cluster(CL).sets);
            for uc = 1:length(unique_clus_subs)
                x = find(STUDY.cluster(CL).sets == unique_clus_subs(uc));
                index = uc;
                if ~isempty(x)
                    if size(x,2)>1 %if subject appears more than once in cluster
                        
                        rm_comp_ind = [rm_comp_ind, x(2:end)];
                    end
                end
            end

new_outlier_comps = STUDY.cluster(CL).comps(rm_comp_ind);
new_outlier_sets =  STUDY.cluster(CL).sets(rm_comp_ind);
cluster_lowestIC(2).comps = [cluster_lowestIC(2).comps, new_outlier_comps];
cluster_lowestIC(2).sets = [cluster_lowestIC(2).sets, new_outlier_sets];
cluster_lowestIC(CL).comps(rm_comp_ind) = [];
cluster_lowestIC(CL).sets(rm_comp_ind) = [];
end
 disp('Subject IC ERSPs consolidated to')
                        disp('one subject/cluster using lowest IC number')
                        
                        
  cd('R:\Ferris-Lab\jacobsen.noelle\Split Belt Pilot Study\FOOF\test')                     
  for CL = clusters_to_plot
      i = i+1;
      all_psd = {psd_results(2).data(CL).data{1, 1}.all_psd}.';
      avg_psd = {psd_results(2).data(CL).data{1, 1}.avg_psd}.';
      all_psd = cat(3, all_psd{:,1});
      all_psd = permute(all_psd,[3 1 2]);
      avg_psd = cat(3, avg_psd{:,1});
      avg_psd = permute(avg_psd,[3 1 2]);
      save(['CL',num2str(CL),'_allpsd'],'all_psd')
      save(['CL',num2str(CL),'_avgpsd'],'avg_psd')
      clear all_psd avg_psd
  end
                        