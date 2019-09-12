function  LegendBoxing_oxDNA1( src,evn,axMain )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% axMain.UserData
       axes(axMain);
       if evn.IntersectionPoint(2)>0
           CreateStruct.Interpreter = 'tex';
           CreateStruct.WindowStyle = 'modal';
           STR={'\fontsize{11}oxDNA interface for generating topology and configuration files: ';
                '(1) Initialize the configuration by clicking the ''initial oxDNA'' button.';
                '(2) Make sure there is no ''?'' for unknown sequences in the caDNAno tab''s table. It may be possible for overhang design but without specifying the closing strands sequences. ';
                '(3) The slider may be helpful for exploding/adjusting the bundles globally.';
                '(4) Fine tuning for rigid-body transformations(RBTs) can be done after exporting configuration and topology files. The operation is similar to Assembly tab.';
                '(5) Specify the indexes of the bundles which are for bended structures. This allows users to break the rigid bundles into three pieces for RBTs to reduce relaxation time.';
                '(6) The results are six files in the directory, including dSRmain.conf for the force file used in relaxation steps only.';
                
           
              
               
               } ;
                                 
                      f = msgbox(STR ,'Instructions', 'help' ,CreateStruct);    

       end

% \color[rgb]{0,0.5,0}


end

