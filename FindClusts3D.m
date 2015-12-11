function C = FindClusts3D (M, D, minClustSz, unifyGap)
% function clusters = FindClusts3D (M, D, minClustSz, unifyGap)
% 
% The function finds the clusters in the 3D boolean matrix 'M'. It starts
% by finding the clusters across the second dimention (usually time), 
% removing small ones (smaller than 'minClustSz') and uniting close ones
% (closer than 'unifyGap'). Then in unites clusteres across the first
% and second dimensions (usually frequency x time). Last, it unites
% clusters according to the adjacency matrix, 'D' across the third
% dimension (usually channels).  
% The function returns a mask on M with 0's where there are no clusters and
% i's where cluster i is.
% 
% Uri Maoz, Caltech. Created 9/17/2014. Major update 10/1/2014

if (~exist('minClustSz','var')||isempty(minClustSz))
    minClustSz=10;
end
if (~exist('unifyGap','var')||isempty(unifyGap))
    unifyGap=5;
end

nCh=size(M,3);
% Find clusters across the first and second dimensions
C=zeros(size(M)); clusts=cell(1,nCh);
for k=1:nCh
    [clusts{k},C(:,:,k)]=FindClusts2D(M(:,:,k),minClustSz,unifyGap);
end
maxChClusts=max(C(:));

iNewClust=1;
nDim=size(C,1)*size(C,2);
if (nCh>1)
    % Connect clusters across the third dimension according to the adjacency
    % matrix 'D'
    iElClusts=find(sum(squeeze(C),1)>0);
    Mconn=cell(nCh,maxChClusts);
    % Run only on electrodes that contain clusters
    for i=iElClusts
        % Connect only clusters of neighboring electrodes
        iNeighbor=find(D(i,:)); iNeighbor=iNeighbor(iNeighbor~=i);
        iNeighbor=iNeighbor(ismember(iNeighbor,iElClusts));
        for j=iNeighbor
            % Intersections between clusters across the 3rd dimension
            [iConn,jConn]=find(C(:,:,i)&C(:,:,j));
            clustFrom=C(sub2ind(size(C),iConn,jConn,repmat(i,size(iConn))));
            clustTo=C(sub2ind(size(C),iConn,jConn,repmat(j,size(iConn))));
            % If clusters have multiple intersections, save only the first
            c=clustFrom+1i*clustTo; c=unique(c);
            clustFrom=real(c); clustTo=imag(c);
            for k=1:length(clustFrom)
                Mconn{i,clustFrom(k)}=[Mconn{i,clustFrom(k)};[j,clustTo(k)]];
            end
        end
    end
    
    % for chFrom=1:size(Mconn,1)
    %  Run only on electrodes with clusters connecting to other electrodes
    floodDone=false;
    if (~isempty(Mconn))
        for chFrom=find(~cellfun(@isempty,Mconn(:,1)))'
            for clustFrom=1:clusts{chFrom}.NumObjects
                if (C((chFrom-1)*nDim+...
                        clusts{chFrom}.PixelIdxList{clustFrom})>0)
                    FloodClust(chFrom,clustFrom,iNewClust);
                    iNewClust=iNewClust+1;
                end
                floodDone=(max(C(:))==0); % True when all clusters were flooded
                if (floodDone), break; end
            end
            if (~isempty(clustFrom)&&floodDone), break; end
        end
    end
end

% if (~floodAny)
%     % if no flooding occurred, then clusters with the same number over
%     % different electrodes (or regions) are actually different clusters
% if (~isempty(Mconn))
%     for iCh=find(~cellfun(@isempty,Mconn(:,1)))'
%         for iClust=1:clusts{iCh}.NumObjects
%             if (C((iCh-1)*nDim+...
%                     clusts{iCh}.PixelIdxList{iClust})>0)
%                 C((iCh-1)*nDim+clusts{iCh}.PixelIdxList{iClust})=-iNewClust;
%             end
%         end
%     end
% end
% At this point, if there are any entries in C that are >0, these are
% clusters that are not connected to any other cluster, and are entirely
% contained within one channel. Find them and give them cluster numbers.
if (any(C(:)>0))
    for iCh=1:size(C,3)
        chC=squeeze(C(:,:,iCh));
        if (any(chC))
            disjointClusts=bwconncomp(chC>0);
            for k=1:disjointClusts.NumObjects
                C((iCh-1)*nDim+disjointClusts.PixelIdxList{k})=-iNewClust;
                iNewClust=iNewClust+1;
            end
        end
    end
%     % We know clusters are disjoint, so they will show up as difference
%     % components on bwconncomp() 
%     disjointClusts=bwconncomp(C>0);
%     for k=1:disjointClusts.NumObjects
%         C(disjointClusts.PixelIdxList{k})=-iNewClust+1-k;
%     end
end
% Turn all negative cluster numbers from flooding to positive
C=-C; 
% % Currently the cluster numbers might be not consecutive, because
% disjoint cluster indices follow all flooded clusters. Make clusters
% consecutively indexed.
u=unique(C(:)); u=u(u~=0);
for k=1:length(u), C(C==u(k))=k; end

% -------------------------------------------------------------------------
    function FloodClust (iCh, iClust, iNewClust)
        % Floodfill from the current cluster to all of the clusters to
        % which it connects
%         if ((iClust<=clusts{iCh}.NumObjects)&&...
%                 (currC(clusts{iCh}.PixelIdxList{iClust}(1))>0))
            C((iCh-1)*nDim+clusts{iCh}.PixelIdxList{iClust})=-iNewClust;
            for iMconn=1:size(Mconn{iCh,iClust},1)
                iNextCh=Mconn{iCh,iClust}(iMconn,1);
                iNextClust=Mconn{iCh,iClust}(iMconn,2);
                if (C((iNextCh-1)*nDim+...
                        clusts{iNextCh}.PixelIdxList{iNextClust})>0)
                    FloodClust(iNextCh,iNextClust,iNewClust);
                end
            end
%         end
    end
% -------------------------------------------------------------------------

end

% ============================= Subfunctions ==============================

function [clusts, C] = FindClusts2D (M, minClustSz, unifyGap)
% function clusts = FindClusts2D (M, minClustSz, unifyGap)
% 
% Find clusters across a 2D matrix 'M', removing clusters less than
% minClustSz large and unifying clusters less than 'unifyGap' apart in one
% of the rows of 'M'. 
% Return the connected-components structure from bwconncomp() and a matrix
% like M with i where cluster i resides
clusts=bwconncomp(M);
% Remove clusters less than minClustSz
iRm=(cellfun(@length,clusts.PixelIdxList)<minClustSz);
clusts.NumObjects=clusts.NumObjects-sum(iRm);
clusts.PixelIdxList(iRm)=[];
clusts.PixelIdxList=... % Transpose all indices lists
    cellfun(@(x) x(:)',clusts.PixelIdxList,'UniformOutput',false);
%Unify clusters less than unifyGap apart
M=false(size(M)); M(cell2mat(clusts.PixelIdxList))=true;
if (clusts.NumObjects>1)
    for i=1:size(M,1)
        if (any(M(i,:)))
            M(i,:)=UnifyClusts1D(M(i,:),unifyGap);
        end
    end
    clusts=bwconncomp(M);
end
C=zeros(size(M)); C=Clusts2Mat(clusts,C);
end

% -------------------------------------------------------------------------
function c = UnifyClusts1D (m, unifyGap)
% Unify clusters along one dimensional array m and return the result in an
% array c with i for cluster i and 0's elsewhere
i1=find(diff(m)>0)+1; if (m(1)), i1=[1,i1]; end
i2=find(diff(m)<0); if (m(end)), i2=[i2,length(m)]; end
iUnif=(i1(2:end)-i2(1:end-1)<=unifyGap);
c=m; 
if (any(iUnif)), c(MultiColon(i2(iUnif),i1([false,iUnif])))=true; end
end

% -------------------------------------------------------------------------
function M = Clusts2Mat (connComp, M)
% The function takes a connected-components structure (output by
% bwconncomp() and adds its components to a 2D matrix, with i for cluster i
for c=1:connComp.NumObjects
    M(connComp.PixelIdxList{c})=c;
end
end

