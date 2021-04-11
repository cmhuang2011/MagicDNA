classdef jsonObject <handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    %  unfinished-      ----CM
    
    
    properties
        filename=[];
        filepath=[];
        Odat=[];
        NumList=[];
        ScafR_AllBase=[];
        Stap_Routing=[];
        StapRoutingCornerRep ;   % C4
        ScafRoutingCornerRep ; % C4
        %----
        OrientBase=10 ;
        period=32 ;
        
        NumListApp=[];
        Stap_Routing_new=[];
        
        Scafdigit=[]
        Stapdigit=[]
        NewDat=[];
        
        %------
        scafSeq=[];
        stapSeq=[];
        stapSegment=[];  %seq
        stapSegment2= []  ; % consider continuity on scaffold
        stapSegR=[];
        ABInds=[];
        ACInds=[];
        TH_seq=[];
        %----
        CylinderConn ;
        
        
        skipPosition=[];
        %------------May 02 2019
        MergedJson_dat=[];
        %------July 22 2020
        DomainStapOnScaf =[] ;
        %----------MagicDNA2 
        skipC4 =[];   %[C4 base] list
        insertC4 =[];   %[C4 base] list
        
        
    end

    
    methods

%         function ExportNewJSON(obj)
%             NNdat=obj.Odat;
%             
% %             prompt = {'Enter File name:'};
% %             dlg_title = 'Input';
% %             num_lines = 1;
% %             defaultans = {'No'};
% %             answer = inputdlg(prompt,dlg_title,num_lines,defaultans);  
%             answer{1} = 'TestJJJ' ;
%             if strcmp(answer,'No')
%                 return;
%             end
%            file_name=strcat(answer{1},'.json');
%            
%                             TTtext=savejson('Title',NNdat2,'ArrayIndent',0,'Compact',1 );
%                 
%                 %                 TTtext=savejson('Title',NNdat,'ArrayIndent',0,'Compact',1 );
%                 TTtext(1:10)=[];
%                 TTtext(end-1:end)=[];
%                 IOfSC2=strfind(TTtext, ',[-999,-888]');  %for color json export
%                 for removedd=1:length(IOfSC2)
%                     UUdataPosittion=strfind(TTtext, ',[-999,-888]');  %for color json export
%                     UUdataPosittion;
%                     Exxtraindex= UUdataPosittion(1);
%                     TTtext(Exxtraindex:Exxtraindex+11)=[];
%                     %                  Cop=strfind(TTtext, ',[-999,-888]');
%                 end
%                 SecTerm999888=strfind(TTtext, '-999,-888');  %for color json export
%                 for removedd=1:length(SecTerm999888)
%                     UUdataPosittion=strfind(TTtext, '-999,-888');  %for color json export
%                     if ~isempty(UUdataPosittion)
%                         Exxtraindex= UUdataPosittion(1);
%                     end
%                     TTtext(Exxtraindex:Exxtraindex+8)=[];
%                 end    
%                 
%                 fileID = fopen([JsonFolder file_name],'w');
%                 fprintf(fileID,TTtext);
%                 fclose(fileID);
%                 fprintf('print json file %s  \n',file_name) ;
%                 
%                 
% %             sdf=3
%             
%         end
        
        
        function obj=jsonObject()
            %         obj=jsonObject()
            %---------debug
            %            [filename,filepath,FilterIndex]= uigetfile({'*.json','JSONfile' },'Select the json file');
            %             dat=loadjson(strcat(filepath,filename));
            %------------------
            
            [obj.filename,obj.filepath,FilterIndex]= uigetfile({'*.json','JSONfile' },'Select the json file');
            dat=loadjson(strcat(obj.filepath,obj.filename));
            obj.Odat=dat;
            
            ColorCode=[-1 ,-1 ,-1];
            EffecCyl=[];
            for k=1:length(dat.vstrands)
                if ~isempty(dat.vstrands{k}.stap_colors)
                    ColorCode=[ColorCode;   [k*ones(size(dat.vstrands{k}.stap_colors,1),1), dat.vstrands{k}.stap_colors   ]];
                    EffecCyl=union(EffecCyl,k);
                end
            end
            ColorCode=setdiff(ColorCode,[-1,-1,-1],'rows') ;
            [~,b]=hist(ColorCode(:,3),unique(ColorCode(:,3)));
            dxy=2.6        ;
            NumList=[-1,-1,-1,-1];    %[k, num, col , row, GcX, GcY]
            Eff=[];
            for k=1:length(dat.vstrands)
%                 if ~ismember(dat.vstrands{k}.num , NumList(:,2) )
                    
                    NumList=[NumList;   [k, dat.vstrands{k}.num,  dat.vstrands{k}.col,  dat.vstrands{k}.row ]    ];
%                 else
%                     MaxN_evv = max(NumList(mod(NumList(:,2),2)==0,2)  ) ;
%                     MaxN_odd = max(NumList(mod(NumList(:,2),2)==1,2)  ) ;
%                     if mod(dat.vstrands{k}.num ,2)==0
%                         NumList=[NumList;   [k, MaxN_evv+2,  dat.vstrands{k}.col,  dat.vstrands{k}.row ]    ];
%                     else
%                         NumList=[NumList;   [k, MaxN_odd+2,  dat.vstrands{k}.col,  dat.vstrands{k}.row ]    ];
%                     end
%                 end
                
            end
            NumList=setdiff(NumList,[-1,-1,-1,-1],'rows');
            NumList(:,end+1:end+2)=0;
            for kk2=1:size(NumList,1)
                NumList(kk2,5)= dxy * NumList(kk2,3) ;
                NumList(kk2,6)= dxy * NumList(kk2,4) ;
            end
            obj.NumList=NumList;
            
            %% get scaf routing
            Scaf_Start=[];  %Cyl-CadIndex
            ColRowStart=[];
            %             ScafR_Skip=ones( 10000,1);
            
            for k=1:length(dat.vstrands)
                ScafInCyl= dat.vstrands{k}.scaf ;
                [~,ind1] =ismember(   ScafInCyl(:,1:2), [-1,-1],'rows') ;
                [~,ind2] =ismember(   ScafInCyl(:,3:4), [-1,-1],'rows') ;
                ind2=~ind2 ;
                All= and (ind1, ind2) ;
                if sum(All)~=0
                    ClyIndx= NumList(NumList(:,1)==k,2) ;
                    Base=find(All)-1;
                    %                 Scaf_Start= [ClyIndx ,Base];
                    Scaf_Start=[Scaf_Start; ones(size(Base))*ClyIndx,Base ];
                    %                 ColRowStart= [dat.vstrands{k}.col,dat.vstrands{k}.row]      ;
                    
                    %---------June 3 2020
                    OnlyColRow = repmat( [ dat.vstrands{k}.col , dat.vstrands{k}.row],  sum(All) ,1) ;
                    ColRowStart=[ColRowStart ;OnlyColRow  ] ;
%                     ColRowStart=[ColRowStart ; dat.vstrands{k}.col , dat.vstrands{k}.row];
                    %                 ScafR_Skip(1)=dat.vstrands{k}.skip(Base+1);
                    %                   ScafR_Skip(1)=0;  %hard
                    %                 break;
                end
            end
            MultiScafR_Skip=ones( 50000,1);   % multi-scaffold , 08/15/2019
            MultiScafAll = zeros(500000,2) ; cScaf =1 ;
            for sf_strandi= 1 :size(Scaf_Start,1)
                ScafR_Skip=ones( 10000,1); % multi
                ScafR_Skip(1)=0;  %hard  % multi
                
                ScafR_AllBase=zeros( 10000,2); ks=2;
                ColRow=zeros( 10000,2);
                ScafR_AllBase(1,:)=Scaf_Start(sf_strandi,:) ;  
                ColRow(1,:)=ColRowStart(sf_strandi,:);
                Current=ScafR_AllBase(1,:) ; nW =0;
                while 1
                    Cyli= NumList(NumList(:,2)==Current(1),1) ;
                    ScafR_AllBase(ks,: ) = dat.vstrands{Cyli}.scaf( Current(2)+1,3:4) ;
                    ks=ks+1 ;
                    ScafR_Skip(ks)= dat.vstrands{Cyli}.skip(Current(2)+1 );
                    ColRow(ks,: ) = [dat.vstrands{Cyli}.col ,dat.vstrands{Cyli}.row] ;
                    Current=dat.vstrands{Cyli}.scaf( Current(2)+1,3:4);
                    if sum(Current==[-1,-1])==2
                        ScafR_AllBase(sum(ScafR_AllBase,2)==0 ,:)=[] ;
                        ColRow(sum(ColRow,2)==0 ,:)=[] ;
                        ScafR_Skip(ScafR_Skip==1)=[];
                        ScafR_AllBase=ScafR_AllBase(1:end-1,:);
                        ColRow=ColRow(2:end,:);
                        ScafR_Skip=ScafR_Skip(2:end,:);
                        break ;
                    end
                    nW=nW+1 ;
                    if  nW==100000
                        fprintf('Error:Circular Scaffold\n')
                        return
                    end
                end
                ScafR_AllBase_Mindex=ScafR_AllBase;
                [~,b2]=ismember(ScafR_AllBase(:,1) , NumList(:,2)) ;
                ScafR_AllBase_Mindex(:,1) = b2;
                obj.ScafR_AllBase{sf_strandi}=ScafR_AllBase ;
                
                MultiScafAll(cScaf:cScaf+size(ScafR_AllBase,1)-1 ,: ) = ScafR_AllBase ;
                MultiScafR_Skip(cScaf:cScaf+size(ScafR_AllBase,1)-1 ,: ) = ScafR_Skip ;
                cScaf=cScaf+size(ScafR_AllBase,1) ;
            end
            MultiScafR_Skip=MultiScafR_Skip(MultiScafAll(:,2)~=0  );  % multi-scaffold , 08/15/2019
            MultiScafAll =MultiScafAll( MultiScafAll(:,2)~=0 ,: ) ;
            skipPosition= MultiScafAll(MultiScafR_Skip==-1,:) ;
            
            %             skipPosition= ScafR_AllBase(ScafR_Skip==-1,:) ;
            
            %             skipPosition(:,2)=skipPosition(:,2)-1 ;
            IndEven= mod(skipPosition(:,1),2)==0 ;
%             skipPosition(IndEven,2)=skipPosition(IndEven,2)-1 ;
%             skipPosition(~IndEven,2)=skipPosition(~IndEven,2)+1 ;
            
            obj.skipPosition=skipPosition;
            
            
            
            %             dfgdg=4
            %% staping routing
            
            m_stpStart=0;
            for k=1:length(dat.vstrands)
                m_stpStart=m_stpStart+   size(dat.vstrands{k}.stap_colors,1);
            end
            Stap_Start=cell(m_stpStart,1 );  %Cyl-CadIndex
            Stap_Routing=cell(m_stpStart,1 );  %Cyl-CadIndex
            
            nC_stap=1 ;
            for k=1:length(dat.vstrands)
                Mat=dat.vstrands{k}.stap_colors ;
                for i=1: size(dat.vstrands{k}.stap_colors )
                    Stap_Start{nC_stap} = [ dat.vstrands{k}.num, dat.vstrands{k}.stap_colors(i,1)] ;
                    nC_stap=nC_stap+1 ;
                    
                end
            end
%             return
            
            for stpi= 1:m_stpStart
                stpi;
                stapI_Rout=zeros( 10000,2); ks=2;
                stapI_Rout(1,:)=Stap_Start{stpi};
                Current=stapI_Rout(1,:) ;
                while 1
                    Cyli= NumList(NumList(:,2)==Current(1),1) ;
                    stapI_Rout(ks,: ) = dat.vstrands{Cyli}.stap( Current(2)+1,3:4) ;
                    ks=ks+1 ;
                    %     ScafR_Skip(ks)= dat.vstrands{Cyli}.skip(Current(2)+1 );
                    %     ColRow(ks,: ) = [dat.vstrands{Cyli}.col ,dat.vstrands{Cyli}.row] ;
                    Current=dat.vstrands{Cyli}.stap( Current(2)+1,3:4) ;
                    if sum(Current==[-1,-1])==2
                        stapI_Rout(sum(stapI_Rout,2)==0 ,:)=[] ;
                        %         ColRow(sum(ColRow,2)==0 ,:)=[] ;
                        %         ScafR_Skip(ScafR_Skip==1)=[];
                        stapI_Rout=stapI_Rout(1:end-1,:);
                        %         ColRow=ColRow(2:end,:);
                        %         ScafR_Skip=ScafR_Skip(2:end,:);
                        break ;
                    end
                end
                
                [isSkip,~]= ismember(stapI_Rout, skipPosition,'rows')  ;
                %             stapI_Rout=stapI_Rout(isSkip==0,:) ;
                Stap_Routing{stpi}=stapI_Rout;
            end
            obj.Stap_Routing=Stap_Routing;
            obj.NumListApp=obj.NumList;
                    obj.getSeq ;
            obj.getSegment;
            obj.getSegment2;
            %--------MagicDNA2, July 30 2020
            skipC4Temp=[-1 -1];
            for k= 1 :length(dat.vstrands)
               Bases = find(dat.vstrands{k}.skip==-1 ) -1 ;  % use python index
               skipC4Temp= [skipC4Temp ; [dat.vstrands{k}.num*ones(length(Bases),1), Bases' ] ] ;
            end 
            skipC4Temp=skipC4Temp(2:end ,:) ;
            obj.skipC4 = skipC4Temp ;
                %----
            insertC4Temp=[-1 -1];
            for k= 1 :length(dat.vstrands)
               Bases = find(dat.vstrands{k}.loop>0 ) -1 ;  % use python index
               insertC4Temp= [insertC4Temp ; [dat.vstrands{k}.num*ones(length(Bases),1), Bases' ] ] ;
            end 
            insertC4Temp=insertC4Temp(2:end ,:) ;
            obj.insertC4 = insertC4Temp ;            
            
        end
        
        function mergeWithOtherJson(obj, otherJSON)
            fprintf('"%s" and "%s" are being merged. \n ',obj.filename ,otherJSON.filename )
            
            JsonName='newMergedSJSON.json';
            OldDat= obj.Odat ;
            NewDat =obj.Odat ;  AddDat= otherJSON.Odat ;
            NewDat.name=JsonName ;
            Zlength1= ceil( size(obj.Odat.vstrands{1}.scaf,1 )/32) *32 ;
            Zlength2= ceil( size(otherJSON.Odat.vstrands{1}.scaf,1 )/32) *32 ;
            
            for k =1 : length(obj.Odat.vstrands)
                ss=size(AddDat.vstrands{k}.stap_colors) ;
                NewDat.vstrands{k}.stap_colors = [ OldDat.vstrands{k}.stap_colors ;  AddDat.vstrands{k}.stap_colors+ Zlength1*ones(ss(1),1)*[1,0] ] ;
                
                tt=size(AddDat.vstrands{k}.stap) ;
                
                AddStap=   AddDat.vstrands{k}.stap+ Zlength1*ones(tt(1),1)*[0,1,0,1]  ;
                Inds=ismember( AddDat.vstrands{k}.stap,[-1,-1,-1,-1] ,'rows') ;
                AddStap(Inds,:) = ones(sum(Inds) ,1)*[-1 -1 -1 -1];
                NewDat.vstrands{k}.stap = [ OldDat.vstrands{k}.stap ; AddStap ] ;
                
                AddScaf=   AddDat.vstrands{k}.scaf+ Zlength1*ones(tt(1),1)*[0,1,0,1]  ;
                Inds2=ismember( AddDat.vstrands{k}.scaf,[-1,-1,-1,-1] ,'rows') ;
                AddScaf(Inds2,:) =ones(sum(Inds2) ,1)*[-1 -1 -1 -1];
                NewDat.vstrands{k}.scaf = [ OldDat.vstrands{k}.scaf ;  AddScaf ] ;
                
                NewDat.vstrands{k}.skip = [ OldDat.vstrands{k}.skip ,  AddDat.vstrands{k}.skip ] ;
                NewDat.vstrands{k}.loop = [ OldDat.vstrands{k}.loop ,  AddDat.vstrands{k}.loop ] ;
                
                
            end
            
            TTtext=savejson('Title',NewDat,'ArrayIndent',0,'Compact',1 );
            
            %                 TTtext=savejson('Title',NNdat,'ArrayIndent',0,'Compact',1 );
            TTtext(1:10)=[];
            TTtext(end-1:end)=[];
            IOfSC2=strfind(TTtext, ',[-999,-888]');  %for color json export
            %                 Cop=IOfSC2;
            %                 for removedd=1:length(IOfSC2)
            %                    Exxtraindex= Cop(1);
            %                      TTtext(Exxtraindex:Exxtraindex+11)=[];
            %                      Cop=strfind(TTtext, ',[-999,-888]');
            %                 end
            
            
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
            
            
            %                 JsonFolder=pwd;
            fileID = fopen([pwd filesep JsonName],'w');
            fprintf(fileID,TTtext);
            fclose(fileID);
            
            
        end
        
        
        function getCylConn(obj)
            XoverBases= find(obj.ScafR_AllBase{1}(1:end-1,1)~=obj.ScafR_AllBase{1}(2:end,1) ) ;
            QQ=obj.ScafR_AllBase{1}(sort([XoverBases;XoverBases+1]),:)  ;
            %            size(QQ) ;
            %            Cyl_s_t= QQ(:,1) ;
            Cyl_s_t=[ QQ(1:2:end,1) ,QQ(2:2:end,1)] +1 ; % avoid 0
            IndSwitch =  Cyl_s_t(:,2)<Cyl_s_t(:,1) ;
            Cyl_s_t(IndSwitch,:) =  flip(  Cyl_s_t(IndSwitch,:) ,2) ;
            %           [EdgeList,b,c]= unique(Cyl_s_t ,'rows') ;
            %           [EdgeCounts,cc2] = hist(c,unique(c)) ;
            
            All_Stap_s_t = [-1,-1] ;
            for k = 1: length(obj.Stap_Routing)
                XoverBases= find(obj.Stap_Routing{k}(1:end-1,1)~=obj.Stap_Routing{k}(2:end,1) ) ;
                QQ=obj.Stap_Routing{k}(sort([XoverBases;XoverBases+1]),:)  ;
                Cyl_s_t_Stap=[ QQ(1:2:end,1) ,QQ(2:2:end,1)] +1 ; % avoid 0
                All_Stap_s_t=[All_Stap_s_t ; Cyl_s_t_Stap ] ;
            end
            All_Stap_s_t(1 , :) =[];
            IndSwitch =  All_Stap_s_t(:,2)<All_Stap_s_t(:,1) ;
            All_Stap_s_t(IndSwitch,:) =  flip(  All_Stap_s_t(IndSwitch,:) ,2) ;
            
            
            [EdgeList,b,c]= unique([Cyl_s_t;All_Stap_s_t] ,'rows') ;
            [EdgeCounts,cc2] = hist(c,unique(c)) ;
            
            for edgei= 1: size(EdgeList,1)
                if  EdgeCounts(edgei) <=2
                    fprintf('Edge btw [%i %i] has Xover by %i times\n' ,EdgeList(edgei,:)-1, EdgeCounts(edgei) )  ;
                end
            end
            
            %           All_Stap_s_t=setdiff(All_Stap_s_t, [-1,-1] ,'rows')
            
            
            
            %           G = graph(Cyl_s_t(:,1),Cyl_s_t(:,2))
            
        end
        
        
        function Out=StapCornerRep(obj)
            CornerRep =cell(size(obj.Stap_Routing)) ;
            for k=1:length(CornerRep)
                BaseRep = obj.Stap_Routing{ k} ;
                Temp=BaseRep(1,:) ;
                for c=2:size(BaseRep,1)
                    if  BaseRep(c,1)==BaseRep(c-1,1) &&  abs(BaseRep(c,2)-BaseRep(c-1,2))==1
                        
                    else
                        Temp= [Temp ; BaseRep(c-1:c,:) ];
                    end
                end
                Temp=[Temp ;   BaseRep(end,:) ];
                CornerRep{k} = Temp ;
            end
            %            sdfs3=3
            obj.StapRoutingCornerRep= CornerRep;
            Out=CornerRep;
        end
        function Out=ScafCornerRep(obj)
            CornerRep =cell(size(obj.ScafR_AllBase)) ;
            for k=1:length(CornerRep)
                BaseRep = obj.ScafR_AllBase{ k} ;
                Temp=BaseRep(1,:) ;
                for c=2:size(BaseRep,1)
                    if  BaseRep(c,1)==BaseRep(c-1,1) &&  abs(BaseRep(c,2)-BaseRep(c-1,2))==1
                        
                    else
                        Temp= [Temp ; BaseRep(c-1:c,:) ];
                    end
                end
                Temp=[Temp ;   BaseRep(end,:) ];
                CornerRep{k} = Temp ;
            end
            %            sdfs3=3
            obj.ScafRoutingCornerRep= CornerRep;
            Out=CornerRep;
        end
        
        %---------------
        function getSeq(obj)
            %           sdfsf=3
            
            
            %           OverHangSeq =   TransFbarOH ;   %case to case
            
            %---------------scaffold domain
            scafSeq=cell(size(obj.ScafR_AllBase)) ;  % 1 by 2
            
            scafOneSeq= strings(size(obj.ScafR_AllBase{1},1  ),1 )  ;  %overhang
            
            
            if length(scafOneSeq)<5
                p8064 =  randseq(length(scafOneSeq))     ;  % random seq
%                p8064 =  p8064Seq   ;   %scaffold sequence   temp to be 8064 
            else
                p8064   =   randseq(length(scafOneSeq))     ;  % random seq
            end
            cc=1;
            for isf= 1 :length(scafOneSeq)
                scafLocation=obj.ScafR_AllBase{1}( isf,:) ;                
                if ismember(scafLocation,  obj.skipPosition , 'rows')
                    scafOneSeq( isf)='' ;
                else
                    scafOneSeq( isf)=p8064(cc) ;
                    cc=cc+1 ;
                end
            end
            obj.scafSeq{1}=scafOneSeq ;
            %----
            if length(obj.ScafR_AllBase) >1
                for k=2: length(obj.ScafR_AllBase)
                    p8064   =   randseq(size(obj.ScafR_AllBase{k},1  ))     ;  % random seq
                    scafTwoSeq= strings(size(obj.ScafR_AllBase{k},1  ),1 )  ;  %8064
                    cc=1;
                    for isf= 1 :length(scafTwoSeq)
                        scafLocation=obj.ScafR_AllBase{k}( isf,:) ;
                        
                        if ismember(scafLocation,  obj.skipPosition , 'rows')
                            scafTwoSeq( isf)='' ;
                        else
                            scafTwoSeq( isf)=p8064(cc) ;
                            cc=cc+1;
                        end
                    end
                    obj.scafSeq{k}=scafTwoSeq;
                end
                ScafALLR= cell2mat((obj.ScafR_AllBase)');
%                 ScafALLSeq=[];
%                 for k=1:length(obj.ScafR_AllBase)
%                 ScafALLSeq=[ScafALLSeq ;obj.scafSeq{1}]   ;
%                 end
%                 ScafALLSeq=[obj.scafSeq{1};obj.scafSeq{2}] ; cell2mat((obj.scafSeq)');
%                 ScafALLSeq=[obj.scafSeq{1};obj.scafSeq{2}] ; 
            else
                ScafALLR=obj.ScafR_AllBase{1};
                ScafALLSeq=scafOneSeq ;
            end
            obj.scafSeq{1}=scafOneSeq;
            
            
            
            ScafALLSeq=  randseq(size(ScafALLR ,1 ))   ;
            
            %---------------------------------
            
            %-----------staple
            obj.Stap_Routing;
            stapSeq=cell(size( obj.Stap_Routing)) ;  % 1 by 2
            
            
            for stpi = 1: length(stapSeq)
                stapR =obj.Stap_Routing{stpi} ;
                strStap=strings(size(stapR,1  ),1 )  ;  %8064
                for BaseOnstran= 1 :size(stapR ,1)
                    
                    Locate= stapR(BaseOnstran ,:) ;
                    [tf, Ind]=ismember(Locate,ScafALLR,'rows');
                    IsSkip =ismember(Locate,  obj.skipPosition , 'rows') ;
                    if tf==1
                        strStap(BaseOnstran)= seqrcomplement(ScafALLSeq(Ind) ) ;
                    elseif IsSkip
                        strStap(BaseOnstran)='';
                    else
                        strStap(BaseOnstran)='?';
                    end
                    
                end
                stapSeq{stpi} = strStap ;
            end
            obj.stapSeq=stapSeq ;
            %             %-------------------
            %             file2_name='testcsvMATLAB.csv ' ;
            %             fileID = fopen( file2_name,'w');
            %             fprintf(fileID , 'Start, End, Sequence, Length, Color  \n' )    ;
            %
            %             for i=1:length(stapSeq)
            %             AA=join(stapSeq{i},'')  ;
            %             fprintf(fileID , '%i[%i], %i[%i], %s, %i,  \n',obj.Stap_Routing{i}(1,:),obj.Stap_Routing{i}(end,:),AA,numel(AA{:})   )    ;
            %             end
            %---------------                   
        end
        
        function getSegment2(obj)
            %            sdfsdf=3
            stapSeg=cell(size(obj.stapSeq)) ;
            stapSegR=cell(size(obj.stapSeq)) ;
            
            ScafR=[];
            for k = 1 : length(obj.ScafR_AllBase)
            ScafR= [ScafR; obj.ScafR_AllBase{k}];
            end
            ScafR =setdiff(ScafR, obj.skipPosition ,'rows','stable' ) ;
            obj.DomainStapOnScaf =cell(size(obj.Stap_Routing ) ) ;
            for is=1:length(obj.Stap_Routing )
                stapi= obj.Stap_Routing{is} ;
                stapi =setdiff(stapi, obj.skipPosition ,'rows','stable' ) ;
                
                [~,IndexOnScaf ]=ismember(stapi , ScafR,'rows') ;
                IndOverhang= find(IndexOnScaf==0) ;
                IndexOnScaf(IndOverhang) = -1:-1:-length(find(IndexOnScaf==0)) ;  % consider overhangs as continuos
                
                MM =  [IndexOnScaf,[0;cumsum(diff(IndexOnScaf)~=-1)]]     ;
                chartable= strings(MM(end)+1 ,1) ;
                StpRcell=cell(MM(end)+1,1) ;
                for segj= 0 :MM(end)
                    Inds=  MM(:,2) ==segj  ;
                    strSS=  join(obj.stapSeq{is}(Inds)','') ;
                    chartable{segj+1}= strSS{1};
                    StpRcell{segj+1}=MM(Inds,:) ;
                end
                stapSeg{is}=chartable ;
                 obj.DomainStapOnScaf{is} =StpRcell ;
                [w,v]=cellfun(@size, StpRcell);
                stapSegR{is} =w' ;
%                                fprintf('Staple %i has segments: %s \n',is,num2str(w')) ;
            end
            %            obj.stapSegment2 =  stapSeg;
            obj.stapSegment2=stapSegR ;
            
        end
        
        
        function getSegment(obj)
            %            sdfsdf=3
            stapSeg=cell(size(obj.stapSeq)) ;
            stapSegR=cell(size(obj.stapSeq)) ;
            for is=1:length(stapSeg )
                MM =  [obj.Stap_Routing{is},[0;cumsum(diff(obj.Stap_Routing{is}(:,1))~=0)]]     ;
                chartable= strings(MM(end)+1 ,1) ;
                StpRcell=cell(MM(end)+1,1) ;
                for segj= 0 :MM(end)
                    Inds=  MM(:,3) ==segj  ;
                    strSS=  join(obj.stapSeq{is}(Inds)','') ;
                    chartable{segj+1}= strSS{1};
                    StpRcell{segj+1}=MM(Inds,:) ;
                end
                stapSeg{is}=chartable ;
                stapSegR{is} =StpRcell ;
            end
            obj.stapSegment =  stapSeg;
            obj.stapSegR=stapSegR ;
            
        end
        
        function visualRouting(obj)
            
            figure(2); clf;  hold on ;
            %------scaffold
            M_scaf=obj.ScafR_AllBase{1}  ; % scaf 1
            [tf,ppY]= ismember(M_scaf(:,1) , obj.NumList(:,2)) ;
            ppY=-ppY ;
            plot(M_scaf(:,2) ,  10*ppY+5 ,'-*b') ;
            
            M_scaf2=obj.ScafR_AllBase{2}  ; % scaf 1
            [tf2,ppY2]= ismember(M_scaf2(:,1) , obj.NumList(:,2)) ;
            ppY2=-ppY2 ;
            plot(M_scaf2(:,2) ,  10*ppY2+5  ,'-.b') ;
            %---------staple
            for stpi= 1 :length( obj.Stap_Routing)
                M_stap= obj.Stap_Routing{stpi}  ; % scaf 1
                [tf,ppY3]= ismember(M_stap(:,1) , obj.NumList(:,2)) ;
                ppY3=-ppY3 ;
                plot(M_stap(:,2) ,  10*ppY3+6 ,'-r') ;
                
            end
            %-----mark skip
            M_skip=obj.skipPosition;
            [tf,yskip]= ismember(M_skip(:,1) , obj.NumList(:,2)) ;
            yskip=-yskip ;
            scatter(M_skip(:,2) ,  10*yskip+5.5 ,'xk') ;
            
        end % end of function visualRouting
        
        function visualStrandMeltT(obj)
            f4H=figure(4); clf;  hold on ;
            %------scaffold
            M_scaf=obj.ScafR_AllBase{1}  ; % scaf 1
            [tf,ppY]= ismember(M_scaf(:,1) , obj.NumList(:,2)) ;
            ppY=-ppY ;
            plot(M_scaf(:,2) ,  10*ppY+5 ,'-*k') ;
            
            M_scaf2=obj.ScafR_AllBase{2}  ; % scaf 1
            [tf2,ppY2]= ismember(M_scaf2(:,1) , obj.NumList(:,2)) ;
            ppY2=-ppY2 ;
            plot(M_scaf2(:,2) ,  10*ppY2+5  ,'-.k') ;
            %---------staple
            minMT=100;  maxmT=0;
            UseMTchoice= 6;   LengthAtleast=8 ;
            for stpi= 1 :length( obj.stapSegR)
                for stpj2 = 1:length(obj.stapSegR{stpi})
                    if length(obj.stapSegment{stpi}{stpj2})> LengthAtleast
                        SeqProperties = oligoprop(obj.stapSegment{stpi}{stpj2} ) ;
                        SeqProperties.Tm(UseMTchoice);
                        obj.stapSegment{stpi}{stpj2};
                        minMT=min([minMT,SeqProperties.Tm(UseMTchoice)]);
                        maxmT=max([maxmT,SeqProperties.Tm(UseMTchoice)]);
                        %                 [minMT, maxmaxmTMT]
                    end
                end
            end
            
            pH2=cell(length( obj.stapSegR) ,1 ) ;
            for stpi= 1 :length( obj.stapSegR)
                M_stap= obj.Stap_Routing{stpi}  ; % scaf 1
                [tf,ppY3]= ismember(M_stap(:,1) , obj.NumList(:,2)) ;
                ppY3=-ppY3 ;
                pH2{stpi}= plot(M_stap(:,2) ,  10*ppY3+6 ,'-r') ;
                MM_MTthissegment=0 ;
                segMT=zeros(length(obj.stapSegR{stpi}),1) ;
                for stpj2 = 1:length(obj.stapSegR{stpi})
                    M_stap= obj.stapSegR{stpi}{stpj2}  ;
                    [tf,ppY3]= ismember(M_stap(:,1) , obj.NumList(:,2)) ;
                    ppY3=-ppY3 ;
                    if length(obj.stapSegment{stpi}{stpj2})> LengthAtleast
                        SeqProperties = oligoprop(obj.stapSegment{stpi}{stpj2} ) ;
                        MTthissegment= SeqProperties.Tm(UseMTchoice) ;
                    else
                        MTthissegment=  minMT;
                    end
                    MM_MTthissegment=max([MM_MTthissegment,     MTthissegment]) ;
                    segMT(stpj2) = MTthissegment ;
                end
                pH2{stpi}.ButtonDownFcn=@(src,evn)showMT_strand(obj,src,evn,f4H ) ;
                pH2{stpi}.UserData.MaxMT=MM_MTthissegment ;
                pH2{stpi}.UserData.segMT=segMT;
                pH2{stpi}.UserData.Ind=stpi ;
            end
            f4H.UserData.pH2=pH2;
            %-----mark skip
            M_skip=obj.skipPosition;
            [tf,yskip]= ismember(M_skip(:,1) , obj.NumList(:,2)) ;
            yskip=-yskip ;
            scatter(M_skip(:,2) ,  10*yskip+5.5 ,'xk') ;
            %--------------
            Str={'--'};
            for k=1:length(pH2)
                Str{k}=strcat('Staple  ',num2str(k));
            end
            popupH = uicontrol('Style', 'popup',...
                'String', Str,'Unit','normalized','Position', [0.8 0.87 0.1 0.08]);
            popupH.Callback=@(src,evn) popupFcn(obj,src,evn ,f4H,pH2   );
            popupH.FontSize=10 ;
            f4H.UserData.popupH=popupH;
            
            
        end % end of function visualStrandMeltT
        
        function popupFcn(obj,src,evn ,f4H  ,pH2 )
            
            showMT_strand(obj,pH2{src.Value},[],f4H  )
        end
        
        
        function showMT_strand(obj,src,evn,f4H  )
            figure(f4H);
            title(  {strcat('Max Melting Temperature = ', num2str(src.UserData.MaxMT,3)  ), strcat('SegMT= ', num2str(src.UserData.segMT',3) )  }    );
            pH2=f4H.UserData.pH2;
            for k=1:length(pH2)
                pH2{k}.Color=[1,0,0] ;
                pH2{k}.LineWidth=0.5;
            end
            set(gca,'FontSize',18);
            
            src.Color= [0,0,1];
            src.LineWidth=1.5;
            f4H.UserData.popupH.Value= src.UserData.Ind;
        end  %end of showMTonTitle
        
        
        
        function checksandwich(obj)
            %         has considered skip
%             Bundle =[];
            for k=1:length( obj.stapSegment)
                
%                 LL=zeros(1, length(obj.stapSegment{k})) ;
%                 for j2= 1:length(obj.stapSegment{k})
%                     LL(j2)= length(obj.stapSegment{k}{j2}) ;
%                 end
%                 LL
                %                fprintf('Staple %i Segment  Lengths= %s \n',k ,num2str(LL) );
                %                stapSegment2
                fprintf('Staple %i, Start %i[%i], L= %i,has segmentL = %s \n',k,obj.StapCornerRep{k}(1,:), sum(obj.stapSegment2{k}) ,num2str(obj.stapSegment2{k})  );
                if sum(obj.stapSegment2{k}) > 60  ||  sum(obj.stapSegment2{k}) < 25
                   fprintf('\n') 
                end
                
%                 if k== 230 
%                     sdsf=3
%                 end
                
            end
            %             sdfsdf=3
        end %end of checksandwich
        
%         function printMT(obj, varargin)
%             %         sdfs=3
%             if nargin>=2
%                 if isnumeric(varargin{1})
%                     stpInds=varargin{1};
%                 elseif strcmp(varargin{1},'All' )
%                     stpInds= 1:length(obj.stapSegment) ;
%                 elseif strcmp(varargin{1},'ModeAB' )
%                     stpInds=[34  35  36  40  41  43  45  47  48  49  50  55  56  57  58  61  65  66  67  68  70  71  72  73  74] ;
%                     obj.ABInds=stpInds ;
%                 elseif strcmp(varargin{1},'ModeAC' )
%                     stpInds=[4    5    8   12   15   23   28   32   33   75   76   77   78   79   80   89   96  105  108  114  123  133  134  135  136  137  140  141  144  145  146  148  152  153  157  158  159] ;
%                     obj.ACInds=stpInds ;
%                 end
%             end
%             %         stpInds
%             %         return
%             
%             UseMTchoice= 6 ; maxT= zeros(length( stpInds),1) ;
%             for i=1:length( stpInds)
%                 k= stpInds(i) ;
%                 TT=zeros(1, length(obj.stapSegment{k})) ;
%                 for j2= 1:length(obj.stapSegment{k})
%                     %             = length(obj.stapSegment{k}{j2}) ;
%                     if length(obj.stapSegment{k}{j2})>8
%                         SeqProperties = oligoprop(obj.stapSegment{k}{j2}) ;
%                         TT(j2)= SeqProperties.Tm(UseMTchoice) ;
%                     else
%                         TT(j2)=0;
%                     end
%                 end
%                 maxT(i)=max(TT) ;
%                 fprintf('Staple %i Segment  MT= %s \n',k ,num2str(TT,3) );
%                 
%             end
%             fprintf('mean max MT = %s \n ', num2str(mean(maxT),3) )
%             fprintf('total print = %i \n ',i )
%             
%             if nargin>=2
%                 if strcmp(varargin{1},'ModeAB' )
%                     figure(22); subplot(2,4,1) ; hold on ;
%                     h1=histogram(maxT,0:1:80,'Normalization','probability' ) ;
%                     title('ModeAB ');
%                     subplot(2,4,2) ;
%                     h1_2=plot(0.5*(h1.BinEdges(1:end-1)+h1.BinEdges(2:end)),cumsum(h1.Values),'b' ) ;
%                     subplot(2,4,[3,4,7,8]) ;hold on ;
%                     h1_2r=plot(0.5*(h1.BinEdges(1:end-1)+h1.BinEdges(2:end)),cumsum(h1.Values),'b' ) ;
%                     
%                     
%                 elseif  strcmp(varargin{1},'ModeAC' )
%                     figure(22); subplot(2,4,5); hold on ;
%                     h3=histogram(maxT,0:1:80,'Normalization','probability'  ) ;
%                     title('ModeAC ');
%                     subplot(2,4,6) ;
%                     h3_2=plot(0.5*(h3.BinEdges(1:end-1)+h3.BinEdges(2:end)),cumsum(h3.Values),'r' ) ;
%                     
%                     subplot(2,4,[3,4,7,8]) ;hold on ;
%                     h3_2r=plot(0.5*(h3.BinEdges(1:end-1)+h3.BinEdges(2:end)),cumsum(h3.Values),'r' ) ;
%                     title(strcat(' MT choice= ', num2str(UseMTchoice)))
%                 end
%             end
%             
%         end %end of printMT
        
        function getTH_seq(obj)
            sacf2=join(obj.scafSeq{1}',''); sacf2=sacf2{1};
            
            TH_seq_AC=cell(length(obj.ACInds ),1) ;
            for k=1:length(obj.ACInds )
                stapInd= obj.ACInds(k) ;
                Fiveprimside=seqrcomplement(obj.stapSegment{stapInd}{1} ) ;
                Threeprimside=seqrcomplement(obj.stapSegment{stapInd}{end}  );
                
                search5p=strfind(sacf2,Fiveprimside) ;
                search3p=strfind(sacf2,Threeprimside) ;
                
                if ~isempty(search5p)
                    if length(search5p)>1
                        take=mod(search5p-1,6)==0 ;
                        search5p=search5p(take);
                    end
                    TH_seq_AC{k}=seqrcomplement(sacf2(search5p:search5p+5));
                    
                elseif  ~isempty(search3p)
                    if length(search3p)>1
                        take=mod(search3p-1,6)==0 ;
                        search3p=search3p(take);
                    end
                    TH_seq_AC{k}=seqrcomplement(sacf2(search3p:search3p+5));
                else
                    for j2=1: length(obj.stapSegment{stapInd})
                        if length(obj.stapSegment{stapInd}{j2})==6
                            Segment=seqrcomplement(obj.stapSegment{stapInd}{j2} ) ;
                            Lca=strfind(sacf2,Segment) ;
                            TH_seq_AC{k}=seqrcomplement(sacf2(Lca:Lca+5));
                            fprintf('staple %i toehold in the middle  \n',stapInd )
                        end
                    end
                    %                    fprintf('toehold on  ' )
                end
                
            end
            obj.TH_seq.TH_seq_AC=TH_seq_AC;
            %-------------------
            TH_seq_AB=cell(length(obj.ABInds ),1) ;
            for k=1:length(obj.ABInds )
                stapInd= obj.ABInds(k) ;
                Fiveprimside=seqrcomplement(obj.stapSegment{stapInd}{1} ) ;
                Threeprimside=seqrcomplement(obj.stapSegment{stapInd}{end}  );
                
                search5p=strfind(sacf2,Fiveprimside) ;
                search3p=strfind(sacf2,Threeprimside) ;
                
                if ~isempty(search5p)
                    if length(search5p)>1
                        take=mod(search5p-1,6)==0 ;
                        search5p=search5p(take);
                    end
                    TH_seq_AB{k}=seqrcomplement(sacf2(search5p:search5p+5));
                    
                elseif  ~isempty(search3p)
                    if length(search3p)>1
                        take=mod(search3p-1,6)==0 ;
                        search3p=search3p(take);
                    end
                    TH_seq_AB{k}=seqrcomplement(sacf2(search3p:search3p+5));
                else
                    s2=0
                end
                
            end
            obj.TH_seq.TH_seq_AB=TH_seq_AB;
        end
        function printTHdeltaG(obj)
            
            dG_AB= zeros(1,length(obj.TH_seq.TH_seq_AB)) ;
            for k=1: length(obj.TH_seq.TH_seq_AB)
                NT=oligoprop(obj.TH_seq.TH_seq_AB{k});
                dG_AB(k)=NT.Thermo(4,3) ;
            end
            %---------
            dG_AC= zeros(1,length(obj.TH_seq.TH_seq_AC)) ;
            for k=1: length(obj.TH_seq.TH_seq_AC)
                NT=oligoprop(obj.TH_seq.TH_seq_AC{k});
                dG_AC(k)=NT.Thermo(4,3) ;
            end
            %-----
            %                   sdf=3
            figure(55);clf;hold on; edges=-10:0.5:-2 ;
            h_AB=histogram(dG_AB,edges,'Normalization','probability');
            h_AC=histogram(dG_AC,edges,'Normalization','probability');
            legend([h_AB,h_AC], 'dG AB','dG AC');
            title(strcat('Mean dG =' , num2str(mean(dG_AB),4),'{  }' ,num2str(mean(dG_AC),4)   ) )
        end
        
        
        function visualMeltT(obj)
            f3H=figure(3); clf;  hold on ;
            %------scaffold
            M_scaf=obj.ScafR_AllBase{1}  ; % scaf 1
            [tf,ppY]= ismember(M_scaf(:,1) , obj.NumList(:,2)) ;
            ppY=-ppY ;
            plot(M_scaf(:,2) ,  10*ppY+5 ,'-*k') ;
            
            M_scaf2=obj.ScafR_AllBase{2}  ; % scaf 1
            [tf2,ppY2]= ismember(M_scaf2(:,1) , obj.NumList(:,2)) ;
            ppY2=-ppY2 ;
            plot(M_scaf2(:,2) ,  10*ppY2+5  ,'-.k') ;
            %---------staple
            minMT=100;  maxmT=0;
            UseMTchoice= 6;   LengthAtleast=8 ;
            for stpi= 1 :length( obj.stapSegR)
                for stpj2 = 1:length(obj.stapSegR{stpi})
                    if length(obj.stapSegment{stpi}{stpj2})> LengthAtleast
                        SeqProperties = oligoprop(obj.stapSegment{stpi}{stpj2} ) ;
                        SeqProperties.Tm(UseMTchoice);
                        obj.stapSegment{stpi}{stpj2};
                        minMT=min([minMT,SeqProperties.Tm(UseMTchoice)]);
                        maxmT=max([maxmT,SeqProperties.Tm(UseMTchoice)]);
                        %                 [minMT, maxmaxmTMT]
                    end
                end
            end
            
            pH=cell(length( obj.stapSegR) ,  length( obj.stapSegR) ) ;
            for stpi= 1 :length( obj.stapSegR)
                for stpj2 = 1:length(obj.stapSegR{stpi})
                    M_stap= obj.stapSegR{stpi}{stpj2}  ;
                    [tf,ppY3]= ismember(M_stap(:,1) , obj.NumList(:,2)) ;
                    ppY3=-ppY3 ;
                    if length(obj.stapSegment{stpi}{stpj2})> LengthAtleast
                        SeqProperties = oligoprop(obj.stapSegment{stpi}{stpj2} ) ;
                        MTthissegment= SeqProperties.Tm(UseMTchoice) ;
                    else
                        MTthissegment=  minMT;
                    end
                    colors= [ maxmT-MTthissegment, 0 , MTthissegment-minMT]/(maxmT-minMT) ;
                    pH{stpi, stpj2}=plot(M_stap(:,2) ,  10*ppY3+6 ,'Color',  colors        ) ;
                    pH{stpi, stpj2}.UserData=MTthissegment ;
                    pH{stpi, stpj2}.ButtonDownFcn=@(src,evn)showMTonTitle(obj,src,evn,f3H ) ;
                    %                 contour(M_stap(:,2) ,  10*ppY3+6,SeqProperties.Tm(end))
                    % sdf=3
                end
            end
            %-----mark skip
            M_skip=obj.skipPosition;
            [tf,yskip]= ismember(M_skip(:,1) , obj.NumList(:,2)) ;
            yskip=-yskip ;
            scatter(M_skip(:,2) ,  10*yskip+5.5 ,'xk') ;
            
        end % end of function visualMeltT
        
        
        
        
        function showMTonTitle(obj,src,evn,f3H)
            figure(f3H);
            title(strcat('Melting Temperature = ', num2str(src.UserData)  ) );
        end  %end of showMTonTitle
        
        
        function exportNewJson(obj)
            
            %-----------
            ArrSize = 384 ;
            
            NewDat= obj.Odat ;
            
            for k = 1:length(NewDat.vstrands)
            NewDat.vstrands{k}.skip= NewDat.vstrands{k}.skip(1:ArrSize) ;
            NewDat.vstrands{k}.loop= NewDat.vstrands{k}.loop(1:ArrSize) ;
            
            NewDat.vstrands{k}.scaf= NewDat.vstrands{k}.scaf(1:ArrSize, :) ;
            NewDat.vstrands{k}.stap= NewDat.vstrands{k}.stap(1:ArrSize ,:) ;
            end
            %-------------
            
            
            TTtext=savejson('Title',NewDat,'ArrayIndent',0,'Compact',1 );
            
            %                 TTtext=savejson('Title',NNdat,'ArrayIndent',0,'Compact',1 );
            TTtext(1:10)=[];
            TTtext(end-1:end)=[];
            IOfSC2=strfind(TTtext, ',[-999,-888]');  %for color json export
            %                 Cop=IOfSC2;
            %                 for removedd=1:length(IOfSC2)
            %                    Exxtraindex= Cop(1);
            %                      TTtext(Exxtraindex:Exxtraindex+11)=[];
            %                      Cop=strfind(TTtext, ',[-999,-888]');
            %                 end
            
            
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
            
            
            JsonFolder=[ pwd filesep];
            fileID = fopen([JsonFolder 'NewJson.json'],'w');
            fprintf(fileID,TTtext);
            fclose(fileID);
            %                 fprintf('print json file %s  \n',file_name) ;
        end
        
        
        %------------------
        %        function findOH(obj)
        %             %%
        %             m_stpStart=length(obj.Stap_Routing) ;
        %             Tol=5; N_Overhanges=0;
        %             TooFarCell=cell(m_stpStart,1) ;
        %             TooFarCell2=cell(m_stpStart,1) ;
        %
        %             for stpi=    1:m_stpStart
        %             CylCenterPos= zeros(size(obj.Stap_Routing{stpi},1) ,3) ; TooFar=CylCenterPos(:,1)==1 ;
        %             [~,CylInd]=ismember(obj.Stap_Routing{stpi}(:,1) , obj.NumList(:,2))  ;
        %             CylCenterPos(:,1:2)= obj.NumList(CylInd , 5:6) ;
        %             CylCenterPos(:,3)= 0.34*obj.Stap_Routing{stpi}(:,2) ;
        %                 for k=1:length(TooFar)-1
        %                     if norm(CylCenterPos(k,:)-CylCenterPos(k+1,:)) >Tol
        %                     TooFar(k)=1 ;
        %                     end
        %                 end
        %             N_Overhanges=N_Overhanges+ sum(TooFar) ;
        %             TooFarCell{stpi}= TooFar ;
        %
        %             TooFarCell2{stpi}=cumsum(circshift(TooFar,1)) ;
        %             end
        %
        % %             [ TooFarCell{4},TooFarCell2{4},obj.Stap_Routing{4}]
        % %             TooFarCell2
        %           ck=1;ck2=1;
        % %  [ TooFarCell{k},TooFarCell2{k},obj.Stap_Routing{k}]
        %             obj.NumListApp=obj.NumList ;
        %             obj.Stap_Routing_new=obj.Stap_Routing ;
        %
        % %             obj
        %              %------changing table
        %             for k=1:m_stpStart
        %                if  sum(TooFarCell2{k}(end-5:end)==1 )==6 && TooFarCell2{k}(end-6) ==0
        %                 [ck,k];
        %                 ck=ck+1;
        %                 CylBase= obj.Stap_Routing{k}(TooFarCell{k}==1,:) ;
        %                 obj=findOHXover(obj, CylBase,k);
        %
        %                elseif sum(TooFarCell2{k}(1:6)==0 )==6 && TooFarCell2{k}(7) ==1
        %                 [ck2,k];
        %                 ck2=ck2+1;
        %                 CylBase= obj.Stap_Routing{k}(find(TooFarCell{k}==1)+1,:) ;
        %                 obj=findOHXover(obj, CylBase,k);
        %                elseif max(TooFarCell2{k})==2
        %                   k ;
        %                end
        %             end
        %             %------changing routing
        %                       ck=1;ck2=1;
        %             for k=1:m_stpStart
        %                if  sum(TooFarCell2{k}(end-5:end)==1 )==6 && TooFarCell2{k}(end-6) ==0
        %                 [ck,k];
        %                 ck=ck+1;
        %                 CylBase= obj.Stap_Routing{k}(TooFarCell{k}==1,:) ;
        %                 obj=findOHXover2(obj, CylBase,k,1);
        %
        %                elseif sum(TooFarCell2{k}(1:6)==0 )==6 && TooFarCell2{k}(7) ==1
        %                 [ck2,k];
        %                 ck2=ck2+1;
        %                 CylBase= obj.Stap_Routing{k}(find(TooFarCell{k}==1)+1,:) ;
        %                 obj=findOHXover2(obj, CylBase,k,2);
        %                end
        %             end
        %
        %
        %
        % %            dfdsf=4
        %             obj
        %
        %        end
        
        %        function obj=findOHXover(obj, CylBase,kthStp)
        %
        %            kRout=CylBase ;
        %            Centers= zeros(size(kRout,1) ,3);
        %            [~,Ind]=ismember(kRout(:,1) , obj.NumList(:,2)) ;
        %            Centers(:,1:2)= obj.NumList(Ind , 5:6) ;
        %            Centers(:,3)= 0.34* kRout(:,2);
        %            QQ=mod(kRout(:,2)-obj.OrientBase,obj.period) + (mod(kRout(:,1),2)==1)*5.5 ;
        %            PhseAng= wrapTo360(QQ/32*3*360) ;
        %            dxy= [cosd(PhseAng) , sind(PhseAng) ] ;
        %
        %            dangle= wrapTo360(atan2d(dxy(2),dxy(1))) ;
        %            RdAng=round(dangle/90)*90 ;
        %            RoundXY=[cosd(RdAng),sind(RdAng)];
        %
        % %            RoundXY= round(dxy);
        %            if kthStp==74
        %                dfsf=3;
        %                RoundXY=[0,-1];
        %            end
        %
        %
        %            Coeffdxy=2.6        ;
        %            NewCylXY= Centers(1:2)+ Coeffdxy*RoundXY ; tol=0.01;
        %            MMNumListApp =obj.NumListApp;
        %            if      ~ismembertol(NewCylXY,obj.NumListApp(:,5:6),tol,'ByRows',true)
        %                if mod(kRout(1),2)==1
        %                    New_R2=  max(MMNumListApp(mod(MMNumListApp(:,2),2)==0,2))+2 ;
        %                else
        %                   New_R2=   max(MMNumListApp(mod(MMNumListApp(:,2),2)==1,2))+2     ;
        %                end
        %           NewCylindLine=[ MMNumListApp(end,1)+1, New_R2 ,RoundXY+ MMNumListApp(Ind,3:4),NewCylXY];
        %           obj.NumListApp=[  MMNumListApp;   NewCylindLine];
        %            else
        %                HaveRepeat=1;
        %            end
        %        end
        %        function obj=findOHXover2(obj, CylBase,StpInd,ccase)
        %
        %            kRout=CylBase ;
        %            Centers= zeros(size(kRout,1) ,3);
        %            [~,Ind]=ismember(kRout(:,1) , obj.NumList(:,2)) ;
        %            Centers(:,1:2)= obj.NumList(Ind , 5:6) ;
        %            Centers(:,3)= 0.34* kRout(:,2);
        %            QQ=mod(kRout(:,2)-obj.OrientBase,obj.period) + (mod(kRout(:,1),2)==1)*5.5 ;
        %            PhseAng= wrapTo360(QQ/32*3*360) ;
        %            dxy= [cosd(PhseAng) , sind(PhseAng) ] ;
        % %            RoundXY= round(dxy);
        %            dangle= wrapTo360(atan2d(dxy(2),dxy(1))) ;
        %            RdAng=round(dangle/90)*90 ;
        %            RoundXY=[cosd(RdAng),sind(RdAng)];
        %            if StpInd==74
        %                dfsf=3;
        %                RoundXY=[0,-1];
        %            end
        %
        %
        %            Coeffdxy=2.6        ;
        %            NewCylXY= Centers(1:2)+ Coeffdxy*RoundXY ; tol=0.01;
        %            MMNumListApp =obj.NumListApp;
        %            if      ~ismembertol(NewCylXY,obj.NumListApp(:,5:6),tol,'ByRows',true)
        %            else
        %          targetXY=RoundXY+ MMNumListApp(Ind,3:4) ;
        %          [~,bq]=ismembertol(targetXY,obj.NumListApp(:,3:4),tol,'ByRows',true)  ;
        %          if ccase==1
        %
        %           obj.Stap_Routing_new{StpInd}(end-5:end,1) =obj.NumListApp(bq,2) ;
        %           AA=obj.Stap_Routing_new{StpInd}(end-6,2);
        %           BB=obj.Stap_Routing_new{StpInd}(end-7,2);
        %
        %           NewZ= linspace(AA,AA+5*(BB-AA),6)';
        %           obj.Stap_Routing_new{StpInd}(end-5:end,2) =NewZ;
        %          elseif ccase==2
        %             sdfsf=3 ;
        %            obj.Stap_Routing_new{StpInd}(1:6,1) =obj.NumListApp(bq,2) ;
        %           AA=obj.Stap_Routing_new{StpInd}(7,2);
        %           BB=obj.Stap_Routing_new{StpInd}(8,2);
        %           NewZ= linspace(AA,AA+5*(BB-AA),6)';
        %           obj.Stap_Routing_new{StpInd}(6:-1:1,2) =NewZ;
        %
        %          end
        %
        %                dosomething=1;
        %            end
        %        end
        
        
        function obj=plotstaple3D(obj)
            figure(3);clf; hold on;
            for k=1:length(obj.Stap_Routing)
                kRout=obj.Stap_Routing{k} ;
                Centers= zeros(size(kRout,1) ,3);
                %            Ind=kRout(:,1);
                [~,Ind]=ismember(kRout(:,1) , obj.NumList(:,2)) ;
                Centers(:,1:2)= obj.NumList(Ind , 5:6) ;
                Centers(:,3)= 0.34* kRout(:,2);
                %            mod(kRout(:,2),obj.period)
                %            dXY=
                obj.OrientBase=10 ;
                QQ=mod(kRout(:,2)-obj.OrientBase,obj.period) + (mod(kRout(:,1),2)==1)*5.5 ;
                PhseAng= wrapTo360(QQ/32*3*360) ;
                dxy= [cosd(PhseAng) , sind(PhseAng) ];
                plot3(Centers(:,1)+dxy(:,1),Centers(:,2)+dxy(:,2),Centers(:,3) ,'-x' )
            end
            xlabel('x') ; ylabel('y') ;zlabel('z') ;
            axis equal
            
        end
        %        function obj=plotstaple3D_new(obj)
        %            figure(3);clf; hold on;
        %            for k=1:length(obj.Stap_Routing_new)
        %            kRout=obj.Stap_Routing_new{k} ;
        %            Centers= zeros(size(kRout,1) ,3);
        % %            Ind=kRout(:,1);
        %            [~,Ind]=ismember(kRout(:,1) , obj.NumListApp(:,2)) ;
        %            Centers(:,1:2)= obj.NumListApp(Ind , 5:6) ;
        %            Centers(:,3)= 0.34* kRout(:,2);
        % %            mod(kRout(:,2),obj.period)
        % %            dXY=
        %             obj.OrientBase=10 ;
        %            QQ=mod(kRout(:,2)-obj.OrientBase,obj.period) + (mod(kRout(:,1),2)==1)*5.5 ;
        %            PhseAng= wrapTo360(QQ/32*3*360) ;
        %            dxy= [cosd(PhseAng) , sind(PhseAng) ];
        %            plot3(Centers(:,1)+dxy(:,1),Centers(:,2)+dxy(:,2),Centers(:,3) ,'-' )
        %            end
        %            xlabel('x') ; ylabel('y') ;zlabel('z') ;
        %            axis equal
        %
        %        end
        
        %        function convertScaf(obj)
        %         MMz=0;
        %         MMz=max([MMz, max(obj.ScafR_AllBase(:,2))]) ;
        %         MMz= 32*(ceil(MMz/32)+1) ;
        %         n_cyl= size(obj.NumListApp,1) ;
        %         SCafDigiM=-1*ones(MMz,4,n_cyl) ;
        %
        %         Current=obj.ScafR_AllBase(1,:) ;
        %         Next=obj.ScafR_AllBase(2,:) ;
        %         IndFirst=find(Current(1)==obj.NumListApp(:,2)) ;
        %         SCafDigiM(Current(2)+1,3:4,IndFirst)=Next ;
        %         Prev=Current ;
        %         for k=2:size(obj.ScafR_AllBase,1)-1
        %         Current=obj.ScafR_AllBase(k,:) ;
        %         Next=obj.ScafR_AllBase(k+1,:) ;
        %         IndFirst=find(Current(1)==obj.NumListApp(:,2)) ;
        %         SCafDigiM(Current(2)+1,3:4,IndFirst)=Next ;
        %         SCafDigiM(Current(2)+1,1:2,IndFirst)=Prev ;
        %           Prev=Current ;
        %         end
        %         SCafDigiM(Next(2)+1,1:2,IndFirst)=Prev;
        %
        %         obj.Scafdigit=SCafDigiM;
        %        end
        
        %        function convertStap(obj)
        %
        %         StapDigiM=-1*ones(size(obj.Scafdigit)) ;
        %         Exclude=[74];
        %         Collect=[-1,-1];
        %         for k= 1:   length(obj.Stap_Routing_new)
        % %             if ismember(k,Exclude)
        % %                 continue
        % %             end
        %
        %         Rout=  obj.Stap_Routing_new{k} ;
        %         Rout_Old=  obj.Stap_Routing{k} ;
        % %         CC=[Rout, Rout_Old];
        % %         QQ=CC(CC(:,1)~=CC(:,3) ,:);
        % %         Collect=[Collect;QQ(:,1:2)];
        % % %         Rout=  obj.Stap_Routing{k} ;
        % %       KK=  mod( unique(Rout(:,1),'stable'),2);
        % %       JJ=KK(2:end)==KK(1:end-1);
        % %       if sum(JJ)~=0
        % %          sdsf=3 ;
        % %       end
        %
        %
        %         Current=Rout(1,:) ;
        %         Next=Rout(2,:) ;
        %         IndFirst=find(Current(1)==obj.NumListApp(:,2)) ;
        %         StapDigiM(Current(2)+1,3:4,IndFirst)=Next ;
        %         Prev=Current ;
        %             for k2=2:size(Rout,1)-1
        %             Current=Rout(k2,:) ;
        %             Next=Rout(k2+1,:) ;
        %             IndFirst=find(Current(1)==obj.NumListApp(:,2)) ;
        %
        %             if sum(StapDigiM(Current(2)+1,1:4,IndFirst))~=-4
        %                 sdfsf=34;
        %                 k
        %             end
        %
        %             StapDigiM(Current(2)+1,3:4,IndFirst)=Next ;
        %             StapDigiM(Current(2)+1,1:2,IndFirst)=Prev ;
        %
        %
        %             Prev=Current ;
        %             end
        %         StapDigiM(Next(2)+1,1:2,IndFirst)=Prev;
        %
        %         end
        %         obj.Stapdigit=StapDigiM;
        %
        %
        % %            Scollect=sortrows(Collect) ;
        % %            figure;plot(Scollect(2:end,2)) ;
        %        end
        
        %        function findStapIndex(obj,CylBaseCadnano,OldOrNew )
        %            if OldOrNew==1
        %                Cells=obj.Stap_Routing  ;
        %            else
        %               Cells= obj.Stap_Routing_new ;
        %            end
        %
        %            for iS=1:length(Cells)
        %               if sum(Cells{iS}(1,:)==CylBaseCadnano)==2
        %                iS
        %                break;
        %               end
        %            end
        %        end
        
        
        
        %        function ExportJSON(obj)
        %             NNdat=loadjson('CadNano.json');
        %             prompt = {'Enter File name:'};
        %             dlg_title = 'Input';
        %             num_lines = 1;
        %             defaultans = {'No'};
        %             answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
        %
        %             if strcmp(answer,'No')
        %             return;
        %             end
        %             file_name=  answer{1};
        %             file_name=strcat(file_name,'.json');%             file_name=JJ;
        %             NNdat.name=file_name;
        %             ColorChoice=[29184      243362     1507550     3355443     5749504     7536862     8947848    11184640    12060012    13369344    16225054];
        %             for cylindex=2:size(obj.Scafdigit,3)
        %                NNdat.vstrands{cylindex} =NNdat.vstrands{1};       %initialized
        %             end
        %             DCellList2=obj.Scafdigit;
        %
        %             ForVisualCylOrder=sortrows(obj.NumListApp,[4,3]) ;
        %             Vec=ForVisualCylOrder(:,1) ;
        %
        %              for cylindex=1: size(obj.Scafdigit,3)
        %                 %                  size(obj.Scafdigit,1)
        %
        %                 NNdat.vstrands{cylindex}.num=obj.NumListApp(Vec(cylindex),2) ;
        %                 NNdat.vstrands{cylindex}.col=obj.NumListApp(Vec(cylindex),3) ;
        %                 NNdat.vstrands{cylindex}.row=obj.NumListApp(Vec(cylindex),4) ;
        %
        %                 ScafM= DCellList2(:,:,Vec(cylindex));
        %                 NNdat.vstrands{cylindex}.scaf=ScafM;
        %
        %                 StapM= obj.Stapdigit(:,:,Vec(cylindex));
        % %                 if cylindex==11
        %                 NNdat.vstrands{cylindex}.stap= StapM ;
        % %                 else
        % %                 NNdat.vstrands{cylindex}.stap= -1*ones(size(ScafM))  ;
        % %                 end
        %                  %----
        %                 NNdat.vstrands{cylindex}.scafLoop  =cell(0,1) ;
        %                 NNdat.vstrands{cylindex}.stapLoop  =cell(0,1)  ;
        %                %----
        %                 SkipVec=zeros(1,size(ScafM,1)) ;
        %                 CylIndCol2= ForVisualCylOrder(cylindex,2) ;
        %                 WithSkipInd=  obj.skipPosition(obj.skipPosition(:,1)==CylIndCol2,2)  ;
        %                 SkipVec(WithSkipInd)=-1;
        %                 NNdat.vstrands{cylindex}.skip  =SkipVec ;
        %                 %----
        %                 NNdat.vstrands{cylindex}.loop  =zeros(1,size(ScafM,1))  ;
        %
        %                 %---
        %
        %                 Head= find(and(StapM(:,3)~=-1,StapM(:,1)==-1))-1 ;
        %                 Qw=ColorChoice( randi(length(ColorChoice),[size(Head,1),1] ) )';
        %                 ColoM=[Head, Qw] ;
        %                 ColoM(end+1,:)=[-999,-888];
        %                 NNdat.vstrands{cylindex}.stap_colors=ColoM ;
        % %                 sdfs=3
        %              end
        %
        %              %------
        %              NNdat2=NNdat ;
        %              TTtext=savejson('Title',NNdat2,'ArrayIndent',0,'Compact',1 );
        %             TTtext(1:10)=[];
        %             TTtext(end-1:end)=[];
        %             IOfSC2=strfind(TTtext, ',[-999,-888]');  %for color json export
        %             for removedd=1:length(IOfSC2)
        %                  UUdataPosittion=strfind(TTtext, ',[-999,-888]');  %for color json export
        %                  UUdataPosittion;
        %                  Exxtraindex= UUdataPosittion(1);
        %                  TTtext(Exxtraindex:Exxtraindex+11)=[];
        % %                  Cop=strfind(TTtext, ',[-999,-888]');
        %             end
        %
        %             SecTerm999888=strfind(TTtext, '-999,-888');  %for color json export
        %             for removedd=1:length(SecTerm999888)
        %                  UUdataPosittion=strfind(TTtext, '-999,-888');  %for color json export
        %                  if ~isempty(UUdataPosittion)
        %                  Exxtraindex= UUdataPosittion(1);
        %                  end
        %                   TTtext(Exxtraindex:Exxtraindex+8)=[];
        %             end
        %
        %            JsonFolder=pwd;
        % %             fileID = fopen(file_name,'w');
        %              fileID = fopen([JsonFolder filesep file_name],'w');
        %             fprintf(fileID,TTtext);
        %             fclose(fileID);
        % %             ReadTest=loadjson(file_name)
        %             ReadTest=loadjson([JsonFolder filesep file_name])
        %              obj.NewDat=ReadTest;
        %              [JsonFolder filesep file_name]
        %
        %        end
        
        
        
    end
    
end

