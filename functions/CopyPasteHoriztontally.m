%this script use to copy and past stap/scaf routing in cadnano "vertically" 

clear ; clc;

[JSON_filename,PathName3,FilterIndex]= uigetfile({'*.json','JSONfile' },'Select the json file');
  dat=loadjson(strcat(PathName3,JSON_filename));
  
  NumList=[-1,-1,-1,-1];    %[k, num, col , row]
   for k=1:length(dat.vstrands)
   NumList=[NumList;   [k, dat.vstrands{k}.num,  dat.vstrands{k}.col,  dat.vstrands{k}.row ]    ];   
  end
 NumList=setdiff(NumList,[-1,-1,-1,-1],'rows');
 
 OldTotalCyl= size(NumList,1) ;
 %--------------------
 %------cylinder index use cadnano, always refers to cadnano;
 
 
%  CopyCyl=[70 ; 77];  %--------------
 
%   CopyCyl=[50;55;52;57;73;70;77   ];  %--------------

%  CopyCyl=[2 ; 3] ; 
%  rangeI =[36, 54 ] ;   

% CC=NumList(1:37 , 2) 
% CC=NumList(25:48 , 2) ;
% CC=NumList(9:16 , 2) 
CC=NumList(34:end-1 , 2) ;

CC=[32,33,36:39]  ;    % the less, the better
CC=NumList(1:end , 2)  ;
%%----------manually
CopyCyl=CC;
%   CopyCyl=[1;0;2;3;5;4;7;6] ; 
 rangeI =[15 , 345 ] ;     %rangeI =[104 , 167 ] ;   

 MoveN_32 =  41 ;   % hard    %should be related to Range2
%------------end of manually
 
 
 Range2=  32*[floor(rangeI(1)/32) ,ceil(rangeI(2)/32)];
% Range2=[0,960] ;
 
 nNew=length(CopyCyl)  ; 
 kCylCopy=zeros(size(CopyCyl)) ;
 for i1=1:length(CopyCyl)
 kCylCopy(i1) =  NumList(NumList(:,2)==CopyCyl(i1),1) ;
 end
 
 
  Ndat =dat ;
   Ndat.name= strcat(  Ndat.name(1:end-5),'_new','.json') ;
   
  for cyli= 1:length( CopyCyl)
   Indk= kCylCopy(cyli) ;
   
   
   stapM= Ndat.vstrands{Indk}.stap ;
%    scafM= Ndat.vstrands{Indk}.scaf ;
   stapSec=  stapM( Range2(1)+1:Range2(2),:) ;
   stapSec2= stapSec;
   stapSec2(:,2)=stapSec2(:,2)+ MoveN_32*32 ; stapSec2( stapSec(:,2) ==-1 ,2) =-1 ;
   stapSec2(:,4)=stapSec2(:,4)+ MoveN_32*32 ;stapSec2( stapSec(:,4) ==-1 ,4) =-1 ;
   stapM( Range2(1)+1+ MoveN_32*32:Range2(2)+ MoveN_32*32,:) = stapSec2 ;
   Ndat.vstrands{Indk}.stap = stapM ;
%-------------------------
   stapColor  =Ndat.vstrands{Indk}.stap_colors ;
   Add=[];
   for k= 1: size(stapColor, 1)
       if stapColor(k,1) >Range2(1) && stapColor(k,1) <Range2(2)
       Add=[Add ; stapColor(k,1)+MoveN_32*32 , stapColor(k,2) ] ;
       end
   end
   stapColor = [ stapColor ; Add ] ;
   stapColor= sortrows(stapColor);
   Ndat.vstrands{Indk}.stap_colors =  stapColor;
   
   %-----------------
      scafM= Ndat.vstrands{Indk}.scaf ;
   scafSec=  scafM( Range2(1)+1:Range2(2),:) ;
   scafSec2= scafSec;
   scafSec2(:,2)=scafSec2(:,2)+ MoveN_32*32 ; scafSec2( scafSec(:,2) ==-1 ,2) =-1 ;
   scafSec2(:,4)=scafSec2(:,4)+ MoveN_32*32 ;scafSec2( scafSec(:,4) ==-1 ,4) =-1 ;
   scafM( Range2(1)+1+ MoveN_32*32:Range2(2)+ MoveN_32*32,:) = scafSec2 ;
   Ndat.vstrands{Indk}.scaf = scafM ;
%---------------
   skipArr =Ndat.vstrands{Indk}.skip ;
   Arrsec= skipArr(Range2(1)+1:Range2(2)) ;
   
   skipArr( Range2(1)+1+ MoveN_32*32:Range2(2)+ MoveN_32*32 ) =Arrsec;
   Ndat.vstrands{Indk}.skip = skipArr;
   
%-----------   
   LoopArr =Ndat.vstrands{Indk}.loop ;
   ArrsecLoop= LoopArr(Range2(1)+1:Range2(2)) ;
   
   LoopArr( Range2(1)+1+ MoveN_32*32:Range2(2)+ MoveN_32*32 ) =ArrsecLoop;
   Ndat.vstrands{Indk}.loop = LoopArr;
 
  end
   
   
   
 %--------------------------------
 
 for k=1:length( Ndat.vstrands)    
     if size( Ndat.vstrands{k}.stap_colors,1)==1
      Ndat.vstrands{k}.stap_colors(end+1,:)= [-999,-888];   
     end
 end
 
 
              TTtext=savejson('Title',Ndat,'ArrayIndent',0,'Compact',1 );
            TTtext(1:10)=[];
            TTtext(end-1:end)=[];
            IOfSC2=strfind(TTtext, ',[-999,-888]');  %for color json export
%             Cop=IOfSC2;
            for removedd=1:length(IOfSC2)
                 UUdataPosittion=strfind(TTtext, ',[-999,-888]');  %for color json export
                 UUdataPosittion;
               Exxtraindex= UUdataPosittion(1);
                 TTtext(Exxtraindex:Exxtraindex+11)=[];
%                  Cop=strfind(TTtext, ',[-999,-888]');
            end
            
            SecTerm999888=strfind(TTtext, '-999,-888');  %for color json export
            for removedd=1:length(SecTerm999888)
                 UUdataPosittion=strfind(TTtext, '-999,-888');  %for color json export
                 if ~isempty(UUdataPosittion)
                 Exxtraindex= UUdataPosittion(1);
                 end
                  TTtext(Exxtraindex:Exxtraindex+8)=[];
            end

            
            
            
             fileID = fopen([PathName3 Ndat.name],'w');
            fprintf(fileID,TTtext);
            fclose(fileID);
%             ReadTest=loadjson(file_name)      
            ReadTest=loadjson([PathName3 Ndat.name])      

 
 
 
 
 
 
 