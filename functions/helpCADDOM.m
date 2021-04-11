function output=helpCADDOM(varargin)
% listing some useful functions which may be used when operating cadDOM
%   Detailed explanation goes here
length(varargin) ;
if isempty(varargin)
    
    fprintf('1. helpCADDOM(1, [''filename''] ) = printSTEPtable : use for print STEP table\n')
    fprintf('2. helpCADDOM(2) = SaveLoadTransMandConn(1); : use for saving TMat and Conn\n')
    fprintf('3. helpCADDOM(3) = SaveAxes : use for copying  the current axes\n')
    fprintf('4. helpCADDOM(4, [''filename'']  ) = print(''-r500'',gcf,''Image4_r500'',''-djpeg'') : use for printing current figure\n')
    fprintf('5. helpCADDOM(5) = trajObject: get/show oxDNA trajectory object and folder \n')
    
    fprintf('6. ss_Assembly= findobj(0,''Tag'',''ss_Assembly'') ; \n')
    fprintf(' GetHyperB= ss_Assembly.UserData.HyperBundle ; \n')
   
    fprintf('7. helpCADDOM(7)=AssignCustomHC : Use this function to assign custom HC crosssection \n')
    fprintf('8. helpCADDOM(8)=AssignCustomSQ : Use this function to assign custom SQ crosssection \n')
    fprintf('9. helpCADDOM(9) : Use this function to save the table data in STEP \n')
    fprintf('10. helpCADDOM(10)=ColorfulRMSF : Use this function to export RMSF as BILD file. Calculate RMSF and Click on figure before the command. \n')
    fprintf('11. helpCADDOM(11) : Change graphics in connectivity graph if have used SaveAxes \n')

    fprintf('12. helpCADDOM(12) = SaveAxes_two : use for copying  the current axes with texts in hidden axes\n')
    fprintf('13. helpCADDOM(13) = ShowScaffssDNA( GetHyperB )  : Show the single strand scaffold. Use for old data missing ssDNA or after cadnano modification. Click on main figure first.\n')
  
    
    fprintf('\n') ;
    fprintf( 'set(gca,''Color'', ''none''); \n') ;
    fprintf('export_fig PrintedFile.png -transparent -r200 \n' )
    
    fprintf('ExportLineToChimera(h) \n' )
    
%     substitute_stap( GetHyperB_S, GetHyperB_T, [16:18], [9:11] )
    fprintf('substitute_stap( GetHyperB_S, GetHyperB_T, [16:18], [9:11] : Bundle indexes ) \n' )
    
    fprintf('getGlobalstapleLength ;  : use when inspecting the staple graph for staple lengths in a group.  \n' )
    fprintf('try defining global minDist_StapleXoverFromTwoSide  for the distance of ignoring staple Xovers for ends.   \n' )
 
   
end

if length(varargin) >= 1
    switch varargin{1}
        case 1
            if length(varargin)==1
                printSTEPtable ;   % give table title and saved file name
            elseif length(varargin)==2
                printSTEPtable(varargin{2}) ;   % give table title and saved file name
            else
                fprintf('Invalid sytax \n')
            end
        case 2
            SaveLoadTransMandConn(1);
        case 3
            SaveAxes ;
        case 4
            if length(varargin)==1
                print('-r500',gcf,'Image4_r500','-djpeg')  ;
            elseif length(varargin)==2
                print('-r500',gcf,num2str(varargin{2}),'-djpeg')  ;
            else
                fprintf('Invalid sytax \n')
            end
        case 5
           if length(varargin)==1
            btn_oxDNATraj1= findobj(gcf,'Tag','btn_oxDNATraj1') ; 
            output=btn_oxDNATraj1.UserData.oxDNA_ex ;
            TrajLocation =btn_oxDNATraj1.UserData.oxDNA_ex.PathName ;         
            fprintf('The trajectory shown here is from %s \n' ,TrajLocation)
           end
        case 7
            AssignCustomHC ;
        case 8
            AssignCustomSQ ; 
        case 9
             t_step=findobj(gcf,'Tag','STEPtable') ;
             STEP_table_Data= t_step.Data;              
             save('STEP_table_Data.mat','STEP_table_Data') ;
        case 10
             ColorfulRMSF ;
        case 11
            h_lines=findobj(gca,'Type','Line') ;
            for k = 1:length(h_lines)
                if  sum(h_lines(k).Color==[1,0,0])==3
                    h_lines(k).LineWidth =2 ;
                end
            end
            
            h_text = findobj(gca,'Type','T') ;
            for k = 1:length(h_text)
                if  sum(h_text(k).Color==[0,0,1])==3  && ~strcmp(h_text(k).String, '0')
                    h_text(k).FontSize = 18 ;
                    uistack(   h_text(k),'top');
                elseif ~strcmp(h_text(k).String, '0')
                    h_text(k).FontSize = 18 ;
                    uistack(   h_text(k),'top');                    
                elseif   sum(h_text(k).Color==[0,0,1])==3
                    
                    h_text(k).Visible = 'off';
                end
            end
        case 12
            SaveAxes_two ;
            
        case 13 
            ss_Assembly= findobj(gcf,'Tag','ss_Assembly') ;
            GetHyperB= ss_Assembly.UserData.HyperBundle ;            
            ShowScaffssDNA( GetHyperB ) ;
            
        case -1
            set(gca,'Color', 'none');    % for making figures only 
            export_fig PrintedFile.png -transparent -r200 
        otherwise
            fprintf('Please enter valid option.!! \n')
    end
end


end

