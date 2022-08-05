clc
clear
close all

N=33;   % number of nodes
K=4;    % number of generations
M=30;   % number of nodes with node

define_constants;
[results, success]=rundcpf('case33_gai');
IJ=results.branch(:,F_BUS:T_BUS);
P_B_IJ=results.branch(:,PF);
R_B_IJ=0;
P_B=zeros(N);  %
for l=1:size(IJ,1)
    if P_B_IJ(l)>=0
        P_B(IJ(l,1),IJ(l,2))=P_B_IJ(l);
    else
        P_B(IJ(l,2),IJ(l,1))=-P_B_IJ(l);
    end
end
P_G=zeros(K,N); % 
for k=1:K
    n=results.gen(k,GEN_BUS);
    P_G(k,n)=results.gen(k,PG);
end
P_L=zeros(M,N); % 
P_L_N=results.bus(:,PD);
m=1;
for j=1:N
    if results.bus(j,PD)~=0
        P_L(m,j)=P_L_N(j);
        m=m+1;
    end
end
for i=1:N
    Gi=find(results.gen(:,GEN_BUS)==i);
    P_Nii(i)=sum(P_B(:,i))+sum(results.gen(Gi,PG));
end
P_N=diag(P_Nii); % 
E_G=[225 525 525 875]';% 
E_N=pinv(P_N-P_B')*P_G'*E_G; % , unit: gCO2/(kWh) 
R_B=P_B'*diag(E_N)/1e3;  % , unit: tCO2/h
for l=1:size(IJ,1)
    if R_B(IJ(l,2),IJ(l,1))>=0
        R_B_IJ(l)=R_B(IJ(l,2),IJ(l,1));
    elseif R_B(IJ(l,1),IJ(l,2))>0
        R_B_IJ(l)=-R_B(IJ(l,1),IJ(l,2));
    end
end
R_L=P_L_N.*E_N/1e3;

q=0;
%% set dirction
for i=1:32
    if results.branch(i,14)<0
        q=results.branch(i,1);
        results.branch(i,1)=results.branch(i,2);
        results.branch(i,2)=q;
    end
end
s=results.branch(:,1);
t=results.branch(:,2);

%开始画图
figure;

G=digraph(s,t); % graph
p=plot(G,'layout','force');

labels= p.NodeLabel;
p.NodeLabel =[];
%location
hold on;
Coord=[
    0,6;   
    0,10;
    0,14;
    0,18;
    0,22;
    0,26;
    0,30;
    0,34;
    0,38;
    0,42;
    0,46;
    0,50;
    0.5,50;
    1,50;
    1,46;
    1,42;
    1,38;
    1,34;
    0.5,10;
    1,10;
    1,14;
    1,18;
    -0.5,14;
    -1,14;
    -1,18;
    -0.5,26;
    -1,26;
    -1,30;
    -1,34;
    -1,38;
    -1,42;
    -1,46;
    -1,50;    
   ];
p.XData=[Coord(:,2)'];
p.YData=[Coord(:,1)'];

%连线

for i=1:13
line([Coord(i,2),Coord(i+1,2)],[Coord(i,1),Coord(i+1,1)],'LineWidth',2,'Color','k');
end
for i=14:17
line([Coord(i,2),Coord(i+1,2)],[Coord(i,1),Coord(i+1,1)],'LineWidth',2,'Color','k');
end
for i=19:21
 line([Coord(i,2),Coord(i+1,2)],[Coord(i,1),Coord(i+1,1)],'LineWidth',2,'Color','k');   
end
for i=23:24
 line([Coord(i,2),Coord(i+1,2)],[Coord(i,1),Coord(i+1,1)],'LineWidth',2,'Color','k');   
end
for i=26:32
 line([Coord(i,2),Coord(i+1,2)],[Coord(i,1),Coord(i+1,1)],'LineWidth',2,'Color','k');   
end
line([10,10],[0,0.5],'LineWidth',2,'Color','k');   
line([14,14],[0,-0.5],'LineWidth',2,'Color','k');   
line([26,26],[0,-0.5],'LineWidth',2,'Color','k'); 
%DGs
line([18,18],[0,-0.25],'LineWidth',2,'Color','b'); 
line([42,42],[0.75,1],'LineWidth',2,'Color','b'); 
line([14,17],[-0.5,-0.5],'LineWidth',2,'Color','b'); 

hold on;
%上色
sz=400;
load D:\CODE\MY\2022-04\E_N.mat;

scatter(p.XData,p.YData,sz,E_N,'filled');

YData=[-0.25;0.75;-0.5];
XData=[18;42;17];
E_N_11=[525;525;875];
scatter(XData,YData,sz,E_N_11,'filled');

h=colorbar('location','east');
 caxis([225 400])
 set(get(h,'label'),'string','En','FontSize',15);
 
 %label
 for i=1:9
    text(p.XData(i)-0.25,p.YData(i),labels(i),'fontsize', 12,'FontName', 'Arial', 'Color','k');
end

 for i=10:length(labels)
    text(p.XData(i)-0.45,p.YData(i),labels(i),'fontsize', 12,'FontName', 'Arial', 'Color','k');
 end
 
for i=1:3
    text(XData(i)-0.55,YData(i),'DG','fontsize', 12,'FontName', 'Arial', 'Color','k');
 end

% legend([c1,c2],'reverse PF','forward PF');

