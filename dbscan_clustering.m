% function labels = dbscan_clustering(data, epsilon, minPts)
%     
%     % DBSCAN Clustering
%     % Inputs:
%     % data - NxM matrix where N is the number of points and M is the number of dimensions
%     % epsilon - the maximum radius of the neighborhood
%     % minPts - minimum number of points required to form a dense region (cluster)
%     % Output:
%     % labels - Nx1 vector of cluster labels (0 indicates noise)
% 
%     n = size(data, 1);  % Number of data points
%     labels = zeros(n, 1);  % Cluster labels initialized to zero
%     clusterIdx = 0;  % Initial cluster index
% 
%     % Loop through each point
%     for i = 1:n
%         if labels(i) ~= 0
%             continue;  % Skip if the point is already labeled
%         end
%         
%         % Find the neighbors of point i
%         neighbors = findNeighbors(i, data, epsilon);
%         
%         if numel(neighbors) < minPts
%             labels(i) = -1;  % Mark as noise if fewer than minPts neighbors
%         else
%             clusterIdx = clusterIdx + 1;  % Increment cluster index
%             expandCluster(i, neighbors, clusterIdx, labels, data, epsilon, minPts);
%         end
%     end
% end
% 
% function neighbors = findNeighbors(idx, data, epsilon)
%     % Find neighbors within epsilon distance of point idx
%     distances = sqrt(sum((data - data(idx, :)).^2, 2));
%     neighbors = find(distances <= epsilon);
% end
% 
% function expandCluster(idx, neighbors, clusterIdx, labels, data, epsilon, minPts)
%     % Expand the cluster starting from point idx
%     labels(idx) = clusterIdx;
%     
%     i = 1;
%     while i <= numel(neighbors)
%         point = neighbors(i);
%         
%         if labels(point) == -1
%             labels(point) = clusterIdx;  % Change noise to border point
%         end
%         
%         if labels(point) == 0
%             labels(point) = clusterIdx;  % Label point as part of the cluster
%             newNeighbors = findNeighbors(point, data, epsilon);
%             
%             if numel(newNeighbors) >= minPts
%                 neighbors = [neighbors; newNeighbors];  % Add new neighbors
%             end
%         end
%         i = i + 1;
%     end
% end

%%

function labels = dbscan_clustering(data, epsilon, minPts)
    n = size(data, 1);  % Number of data points
    tempLabels = cell(n, 1);  % Temporary cell array to store labels
    clusterIdx = 0;  % Initial cluster index

    % Set up a parallel pool (if not already active)
    if isempty(gcp('nocreate'))
        parpool;  % Start parallel pool
    end
    
    parfor i = 1:n
        % Initialize labels for this iteration
        tempLabel = zeros(n, 1);

        if any(tempLabel(i) ~= 0)
            tempLabels{i} = tempLabel;
            continue;  % Skip if the point is already labeled
        end
        
        % Find the neighbors of point i
        neighbors = findNeighbors(i, data, epsilon);
        
        if numel(neighbors) < minPts
            tempLabel(i) = -1;  % Mark as noise if fewer than minPts neighbors
        else
            % Increment cluster index (needs to be handled differently)
            tempClusterIdx = clusterIdx + 1;  % Use temporary cluster index
            expandCluster(i, neighbors, tempClusterIdx, tempLabel, data, epsilon, minPts);
        end

        tempLabels{i} = tempLabel;  % Store the result
    end
    
    % Combine results from all iterations
    labels = zeros(n, 1);
    for i = 1:n
        labels = max(labels, tempLabels{i});
    end
    
    % Optionally shut down the parallel pool
    % delete(gcp('nocreate'));
end

function neighbors = findNeighbors(idx, data, epsilon)
    % Find neighbors within epsilon distance of point idx
    distances = sqrt(sum((data - data(idx, :)).^2, 2));
    neighbors = find(distances <= epsilon);
end

function expandCluster(idx, neighbors, clusterIdx, labels, data, epsilon, minPts)
    % Expand the cluster starting from point idx
    labels(idx) = clusterIdx;
    
    i = 1;
    while i <= numel(neighbors)
        point = neighbors(i);
        
        if labels(point) == -1
            labels(point) = clusterIdx;  % Change noise to border point
        end
        
        if labels(point) == 0
            labels(point) = clusterIdx;  % Label point as part of the cluster
            newNeighbors = findNeighbors(point, data, epsilon);
            
            if numel(newNeighbors) >= minPts
                neighbors = [neighbors; newNeighbors];  % Add new neighbors
            end
        end
        i = i + 1;
    end
end
