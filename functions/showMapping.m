
function showMapping(src,evn,StapToScaf_Corr_NBase,axx ,ax ,Q ,CircularH , varargin )
% src

for k=1: length(CircularH)
    CircularH{k}.LineWidth =0.2 ;
    CircularH{k}.Color =[CircularH{k}.Color(1:3), 0.1 ];
end
SelectedStap=[];
axes(axx) ;
AxCP = get(gca,'CurrentPoint') ; 
%  size(StapToScaf_Corr_NBase)
QerryScafStap= round(AxCP(2,1:2) ) ;

if QerryScafStap(1)> 0 && QerryScafStap(1)<= size(StapToScaf_Corr_NBase,2)  && QerryScafStap(2)>= 0 && QerryScafStap(2)<= size(StapToScaf_Corr_NBase,1)
    Value =  StapToScaf_Corr_NBase(QerryScafStap(2) , QerryScafStap(1)) ;
    title(strcat('Scaf =', num2str(QerryScafStap(1)),' Stap = ',   num2str(QerryScafStap(2)), ' Comp. Base =', num2str(Value) )) ;
    CircularH{round(QerryScafStap(2))}.LineWidth =3 ;
    CircularH{round(QerryScafStap(2))}.Color =[CircularH{QerryScafStap(2)}.Color(1:3), 1 ];
      
    axes(ax) ;
    title(strcat(' Stap# ', ' ' , num2str(QerryScafStap(2)) , ' connects to ' ,' ' ,num2str(Q(QerryScafStap(2))), ' scaf(s).'   ) ) ;
    xlabel( sprintf('Mapping [%s]',num2str(StapToScaf_Corr_NBase(QerryScafStap(2),: ) )) ) ;
    SelectedStap =QerryScafStap(2) ;
end

axes(ax) ;
AxCP2 = get(gca,'CurrentPoint')  ;
QerryScafStap= ceil(AxCP2(2,1:2) )  ;
if QerryScafStap(1)> 0 && QerryScafStap(1)<=  1  && QerryScafStap(2)> 0 && QerryScafStap(2)<= size(StapToScaf_Corr_NBase,1)
    title(strcat(' Stap# ', ' ' , num2str(QerryScafStap(2)) , ' connects to ' ,' ' ,num2str(Q(QerryScafStap(2))), ' scaf(s).'   ) ) ;
    %     Value =  StapToScaf_Corr_NBase(QerryScafStap(2) , QerryScafStap(1)) ;
    %     title(strcat('Scaf =', num2str(QerryScafStap(1)),' Satp = ',   num2str(QerryScafStap(2)), ' Comp. Base =', num2str(Value) )) ;
    %     CircularH{round(QerryScafStap(2))}.LineWidth =10 ;
    CircularH{round(QerryScafStap(2))}.LineWidth =3 ;
    CircularH{round(QerryScafStap(2))}.Color =[CircularH{QerryScafStap(2)}.Color(1:3), 1 ];
    
    xlabel( sprintf('Mapping [%s]',num2str(StapToScaf_Corr_NBase(QerryScafStap(2),: ) )) ) ;
     SelectedStap =QerryScafStap(2) ;
else
%         notTrigger=1
   
end
% varargin
% nargin
if ~isempty(varargin{1})
%        varargin : t_json, pStapleH,  plotH ,jsonSlider2
       t_json=varargin{1}{1} ;  plotH=varargin{1}{2} ;  pStapleH=varargin{1}{3} ;  jsonSlider2=varargin{1}{4} ;
    if ~isempty(SelectedStap)
        for k=1:length(pStapleH)
            if ismember(k,SelectedStap )
                pStapleH{k}.LineWidth =5;        plotH{k}.LineWidth =5 ;
                pStapleH{k}.Color(4) = 1 ;      plotH{k}.Color(4) = 1 ;
            else
                pStapleH{k}.LineWidth =2 ;       plotH{k}.LineWidth =2 ;
                pStapleH{k}.Color(4) =jsonSlider2.Value ;  plotH{k}.Color(4)=jsonSlider2.Value ;
            end
        end
    end
end
drawnow;

end