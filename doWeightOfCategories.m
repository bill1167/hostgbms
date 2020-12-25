function [ weightOfCategories, STop, WTop, ITop, CTop ] = doWeightOfCategories( disMat, totalSeq, categories, numberOfClusters, clusterNames, mydisMatCategoryFile, mode)
%%%%%%%%%%%%%%%%%%%% self added, top 5 and their corresponding categories
if (nargin<7)
    mode = 0;
end

S = {};
I = {};
for i=1:totalSeq
    [Stemp{i}, Itemp{i}] = sort(disMat(:,i)); %disMat(1:10,1:10) is a smaller matrix
    S = [S Stemp{i}];
    I = [I Itemp{i}];
end
% top 5 similarity except itself
top = 6;
STop = {};
WTop = {};
ITop = {};
CTop = {};
for i=1:totalSeq
    for j=1:top
        STopTempInfo = S{i}(j,1);
        STopTemp{j} = STopTempInfo;
        if j > 1 % compute weight except itself
            WTopTemp{j} = 1/STopTemp{j};
        else
            WTopTemp{j} = 0;
        end
        ITopTempInfo = I{i}(j,1);
        ITopTemp{j} = ITopTempInfo;
        CTopTempInfo = categories(I{i}(j,1));
        CTopTemp{j} = CTopTempInfo{1};
    end
    STop{i} = cell2mat(STopTemp); %score top 6
    WTop{i} = cell2mat(WTopTemp); %weight top 6
    ITop{i} = cell2mat(ITopTemp); %index top 6
    CTop{i} = CTopTemp; %category top 6
end
% weight of categories
weightOfCategories = {};

deltaWTop = {};
for i=1:totalSeq
    for j=2:top % exclude itself
        deltaWTop{i}(1,j) = WTop{i}(1,j) - WTop{i}(1,top);
    end
    sumDeltaWTop = sum(deltaWTop{i}(1,2:top));
    weightOfRowCategories = [0, 0, 0, 0, 0, 0, 0, 0, 0]; %cell(1, numberOfClusters);
    for j=2:top % exclude itself
        for n=1:numberOfClusters
            if string(CTop{i}(1,j)) == string(clusterNames(n))
                weightOfRowCategories(n) = weightOfRowCategories(n) + deltaWTop{i}(1,j);
            end
        end
    end
    weightOfRowCategories = weightOfRowCategories/sumDeltaWTop;
    weightOfCategories{i} = weightOfRowCategories;  
end

% for i=1:totalSeq
%     rowSum = sum(WTop{i}(1,2:top));
%     weightOfRowCategories = [0, 0, 0, 0, 0, 0, 0, 0, 0]; %cell(1, numberOfClusters);  
%     for j=2:top % exclude itself
%         for n=1:numberOfClusters
%             if string(CTop{i}(1,j)) == string(clusterNames(n))
%                 weightOfRowCategories(n) = weightOfRowCategories(n) + WTop{i}(1,j);
%             end
%         end
%     end
%     weightOfRowCategories = weightOfRowCategories/rowSum;
%     weightOfCategories{i} = weightOfRowCategories;    
% end

% reshape from [1 437] to [437 1] before output. 
% %access numbers are sorted by A1, A10, A11, ..., A19, A2, A20, ...in reshape(AcNmb, [437 1])
if mode == 0 %overwrite existing file
    writetable(cell2table(reshape(weightOfCategories,[totalSeq 1])), mydisMatCategoryFile)
else %add content to existing file
    writetable(cell2table(reshape(weightOfCategories,[totalSeq 1])), mydisMatCategoryFile, 'WriteMode','Append')
end
%%%%%%%%%%%%%%%%%%%%
