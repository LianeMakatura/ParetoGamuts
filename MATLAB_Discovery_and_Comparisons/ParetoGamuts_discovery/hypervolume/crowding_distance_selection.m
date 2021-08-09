function [cd_select, cd_discard] = crowding_distance_selection(F, mu)

	% Select based on non-dominated sorting
	[nds_ranks] = non_dominated_sorting(F);
	i = 1;
	num_selected = 0;
	selected_flags = zeros(1,size(F,2));
	while num_selected < mu
		critical_layer_indexes = find(nds_ranks == i);
		l = length(critical_layer_indexes);
		if (l > 0 && num_selected + l <= mu)
			cd_select(num_selected+1:num_selected+l) = critical_layer_indexes;
			selected_flags(critical_layer_indexes) = 1;
		end
		num_selected = num_selected + l;
		i = i + 1;
	end

	% Select based on crowding distance
	if (num_selected > mu)
		[crowding_distances] = compute_crowding_distances(F(:,critical_layer_indexes));
		[cdsort, cdindex] = sort(crowding_distances, 'descend');
		cd_select((num_selected-l)+1:mu) = critical_layer_indexes(cdindex(1:mu-(num_selected-l)));
		selected_flags(critical_layer_indexes(cdindex(1:mu-(num_selected-l)))) = 1;
	end
	cd_discard = find(selected_flags == 0);
	
end
