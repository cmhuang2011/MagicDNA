function ChangeAxesLimit(src,evn)

% %     set(gcf,'KeyPressFcn',@(src,evn)ChangeAxesLimit(src,evn) )

ax=gca; fH=gcf;
interval =5; R=1.1;
switch  evn.Character
    
    case 'h'
        switch fH.UserData.MovieStep
            case 2
                implay('EditSketchSTEP.mp4');
            case 3
                implay('LineToBundle.mp4');
            case 4
                opts.Interpreter = 'tex'; opts.Default = 'Scaffold algorithm and GUI.';
                 answer = questdlg('\fontsize{15} Which movie do you want to play? ?'  , ...
                'Inspect or not?', ...
                'Scaffold algorithm and GUI.','The options for scaffold algorithm.',opts);
                switch answer
                    case 'Scaffold algorithm and GUI.'
                        implay('Scaffold.mp4');
                    case 'The options for scaffold algorithm.'
                        implay('ScafOption_LowRes.mp4');
                end
            case 5
                implay('Staple.mp4');
            case 6
                implay('StapleOverhang.mp4');
            case 7
                implay('oxDNA1.mp4');
            case 8
                implay('oxDNA2.mp4');
                
        end
    case 'G'
        CalculateOxDNAAngle ;
    case 'g'
        CalculateOxDNAAngle;
        
        
    case 'q'
        ax.XLim= ax.XLim + interval ;
    case 'a'
        ax.XLim= ax.XLim - interval ;
    case 'Q'
        ax.XLim= R*(ax.XLim -mean(ax.XLim))+ mean(ax.XLim)   ;
    case 'A'
        ax.XLim= (ax.XLim -mean(ax.XLim))/R + mean(ax.XLim)   ;
        %-----------
    case 'w'
        ax.YLim= ax.YLim + interval ;
    case 's'
        ax.YLim= ax.YLim - interval ;
    case 'W'
        ax.YLim= R*(ax.YLim -mean(ax.YLim))+ mean(ax.YLim)   ;
    case 'S'
        ax.YLim= (ax.YLim -mean(ax.YLim))/R + mean(ax.YLim)   ;
        %------------
    case 'e'
        ax.ZLim= ax.ZLim + interval ;
    case 'd'
        ax.ZLim= ax.ZLim - interval ;
    case 'E'
        ax.ZLim= R*(ax.ZLim -mean(ax.ZLim))+ mean(ax.ZLim)   ;
    case 'D'
        ax.ZLim= (ax.ZLim -mean(ax.ZLim))/R + mean(ax.ZLim)   ;
        
        %---------------
    case 'x'
        axis equal;
        fprintf('set axis equal \n')
    case 'X'
        axis auto
        fprintf('set axis auto \n')
    case 'P'
        fprintf('Capturing..... \n') ;
        print('-r100',gcf,'SnapShot_r100','-dpng')  ;
        fprintf('Captured current figure. \n') ;
    case 'p'
        fprintf('Capturing..... \n') ;        
        print('-r100',gcf,'SnapShot_r100','-dpng')  ;
        fprintf('Captured current figure. \n') ;
    case 'n'
        axis normal
        fprintf('set axis normal \n')
    case 'N'
        axis normal
        fprintf('set axis normal \n')
%     case 'T'
%         fprintf('doing something .....\n')
%         ss_Assembly= findobj(gcf,'Tag','ss_Assembly') ;
%         GetHyperB= ss_Assembly.UserData.HyperBundle ;
%         if length(GetHyperB.ScafRouting) ==1   % start with one cycle case
%             NOriGPairBCB = GetHyperB.ScafRouting;
%             BaseRoutOri = interpolateBase_ThreeCol( NOriGPairBCB{1} ); NBaseOri = size(BaseRoutOri,1) ;
%           
%             [OneXover, XoverList]=Given2CylceFindXoverList(GetHyperB,NOriGPairBCB{1},NOriGPairBCB{1},[]) ;
%               fprintf('check applying %i Xovers to split the scaf routing and match the range ( 5 percent)   .....\n', size(XoverList,1))
%               fprintf('NBaseOri~= %i    .\n' ,NBaseOri ) ;
%             for k = 1:1 :  size(XoverList,1)
%                 OneXoverFromList= [XoverList(k,1:6) ; XoverList(k,7:12)] ;
%                 TwoCycles=GetHyperB.AddXover2OneCycle(NOriGPairBCB,OneXoverFromList) ;
%                 
%                 BaseRout1 = interpolateBase_ThreeCol( TwoCycles{1} ); NBase1 = size(BaseRout1,1) ;
%                 BaseRout2 = interpolateBase_ThreeCol( TwoCycles{2} ); NBase2 = size(BaseRout2,1) ;
%                 
%                 if NBase1>NBaseOri*0.48 && NBase1<NBaseOri*0.52 &&  NBase2>NBaseOri*0.48 && NBase2<NBaseOri*0.52
%                      fprintf('NBase1~= %i ,  NBase2~= %i  .\n' ,NBase1,NBase2 )
%                     GetHyperB.ScafRouting=TwoCycles ;
%                     ConvertScafG(GetHyperB);
%                     ConvertScafSQ(GetHyperB);      % Get properties:ScafdigitSQ
%                     ConvertScafHC(GetHyperB);      % Get properties:ScafdigitHC
%                     fprintf('Found the crossover satisfying the ranges, update data .....\n') ;
%                     k
%                     break ;
%                 end
%             end
%         else
%                     fprintf('The number of loops is not one. Nothing happens. \n')
%         end

end