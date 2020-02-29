% Simple Random Walk for Fracture Approximation
% Developed by Arash Rabbani rabarash@yahoo.com

addpath('Functions');
clear; clc; close all
S=[250,250]; % Size of the Geometry in pixels
Thickness=3; % Thickness of the fracture pixels
StartLocation=[round(S(1)*.5),1]; % middle of the sample, left side

% Select the input image as the cost map
A=imread('Data/A1.png');
% A=imread('Data/A2.png');
% A=imread('Data/A3.png');
% A=imread('Data/A4.png');


A=imresize(A,S); if ndims(A)==3; A=rgb2gray(A); end;  A=mat2gray(A); A=flatout(A);
CostMap=normal(A);  
Population=1;
ShowSteps=50;
ConvergeCriteria=30;
S=size(CostMap);
Generations=1000;  % Total Number of Generations
[X,Y] = meshgrid([-1:1],[-1:1]); Movs=[X(:) Y(:)]; Dirs=size(Movs,1); % number of possible directions

MinCost=zeros(Generations,1);
Loc=.5;

% This matrix shows the bias of the walkers to move in a specific direction
% If you do not put a large value for right direction, walkers may not
% reach the end of the cost map. We have assumed that they start from the
% left side of geometry and move toward the right side

% By changing top-right and bottom right values in Bias matrix you are able
% to make inclined fractures 

Bias=[0   1   1;... %   top-left    top     top-right
      0   0   2;... %   left       center   right
      0   1   1]';  %   bot-left   bottom   bot-right
Bias=Bias(:);
MaxMoves=S(2)*3; % Maximum number of steps
Best=1e10;
figure;      set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
for I=1:Generations
    Cost=0;
    MoveList=randw(Bias,MaxMoves)'; % randw function generates weighted random numbers that shows different directions 
    if I>1
%     MoveList=(BestMoveList.*.95+MoveList.*.05); MoveList=round(MoveList); MoveList(MoveList==0)=1;
    end
    Exit=0;
    % Starting point
    L=StartLocation;
    Path=zeros(1,MaxMoves); % Pathways in each generation
    id=sub2ind(S,L(:,1),L(:,2));
    Path(:,1)=id;
    for J=2:MaxMoves
        L2=L+Movs(MoveList(:,J),:);
        L2(L2(:,2)<=0,2)=L(L2(:,2)<=0,2); % Freeze on left border
        L2(L2(:,2)>S(2),2)=S(2); % Freeze on right border
        L2(L2(:,1)>S(1),1)=L(L2(:,1)>S(1),1);  % Avoid bottom border
        L2(L2(:,1)<=0,1)=L(L2(:,1)<=0,1);  % Avoid top border
        L(logical(1-Exit),:)=L2(logical(1-Exit),:);
        id=sub2ind(S,L(:,1),L(:,2));
        Exit=L(:,2)>=S(2);
        Cost=Cost+CostMap(id).*(1-Exit);
        Path(:,J)=id;
    end
    Cost=Cost./S(2); 
%     Cost(logical(1-Exit))=mean(CostMap(:)).*2;
    [~,p]=sort(Cost);
    if Best>Cost
        BestMoveList=MoveList;
        Best=Cost;
        B=Path(p(1),:);
        A=CostMap.*0; A(B)=1;
        Mask=1-imdilate(A,ones(Thickness));
        subplot(1,2,1); imagesc(CostMap.*Mask); axis equal tight; c=colorbar; c.Label.String = 'Cost Values';
    end
    BestList(I)=Best;
        if I>1; subplot(1,2,2); hold on; plot([I-1 I],BestList(I-1:I),'k'); set(gca, 'XScale', 'log'); xlabel('Generations'); ylabel('Cost'); drawnow; end
end
imwrite(CostMap.*Mask,'Fractured.png');
