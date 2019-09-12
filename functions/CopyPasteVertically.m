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
 %%
 
%  CopyCyl=[70 , 77];  %--------------
 
%   CopyCyl=[50;55;52;57;73;70;77   ];  %--------------

  CC=  NumList( 1:end ,2) ;
%   CC(CC==43)=[];
  CopyCyl=CC;
 
 
 nNew=length(CopyCyl)  ; 
 kCylCopy=zeros(size(CopyCyl)) ;
 for i1=1:length(CopyCyl)
 kCylCopy(i1) =  NumList(NumList(:,2)==CopyCyl(i1),1) ;
 end
 
 NewCylInd=zeros(length(CopyCyl),4  );
 
 ColShift=0;
 RowShift= 12 ;
 for k=1:length(CopyCyl)
    NewCylInd(k,1)=OldTotalCyl+k ;
    if mod(CopyCyl(k),2)==0 %even
        OldEven= NumList(mod(NumList(:,2),2)==0,2) ;
        OldEven=union(OldEven ,NewCylInd(:,2) ) ;
        OldEven( mod(OldEven,2)==1)=[];
         NewCylInd(k,2) = max(OldEven)+2 ;
    else
        OldEven= NumList(mod(NumList(:,2),2)==1,2) ;
        OldEven=union(OldEven ,NewCylInd(:,2) ) ;
        OldEven( mod(OldEven,2)==0)=[];

         NewCylInd(k,2) = max(OldEven)+2 ;
    end
    
    
     NewCylInd(k,3) =NumList( NumList(:,2)==CopyCyl(k) ,3)+ColShift;
      NewCylInd(k,4) =NumList( NumList(:,2)==CopyCyl(k) ,4)+RowShift ;
   
      if ismember(NewCylInd(k,3:4) ,  NumList(:,3:4) , 'rows')
         showerror=1
          
      end
      
 end
 
 Ndat =dat ;
  Ndat.name= strcat(  Ndat.name(1:end-5),'_CPV','.json') ;
 
 for jj2= 1: length(CopyCyl)
 
  Ndat.vstrands{ NewCylInd(jj2,1)} =   dat.vstrands{ kCylCopy(jj2)};
  
  Ndat.vstrands{ NewCylInd(jj2,1)}.col =NewCylInd(jj2,3) ;
    Ndat.vstrands{ NewCylInd(jj2,1)}.row =NewCylInd(jj2,4) ;

 end
 NewCylInd; 
 
 substable=[CopyCyl , NewCylInd(:,2) ] ;
 
 for upd= 1: size(NewCylInd,1) 
    Indk =  NewCylInd(upd,1) ;
     
   stapM= Ndat.vstrands{Indk}.stap ;
   scafM= Ndat.vstrands{Indk}.scaf ;
   
   for subs= 1: size(substable,1)
     stapM( stapM(:,1) ==substable(subs,1) ,1) =substable(subs,2) ;
     stapM( stapM(:,3) ==substable(subs,1) ,3) =substable(subs,2) ;
     
     scafM( scafM(:,1) ==substable(subs,1) ,1) =substable(subs,2) ;
     scafM( scafM(:,3) ==substable(subs,1) ,3) =substable(subs,2) ;
    
   end
   Ndat.vstrands{Indk}.stap  =   stapM ;
   Ndat.vstrands{Indk}.scaf =  scafM ;
   
   Ndat.vstrands{Indk}.num=  NewCylInd(upd,2) ;
 end
 
 
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

 
 
 
 
 
 
 