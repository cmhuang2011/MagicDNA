function  LegendBoxing_customHC( src,evn,axMain )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% axMain.UserData
       axes(axMain);
       if evn.IntersectionPoint(2)>0
           CreateStruct.Interpreter = 'tex';
           CreateStruct.WindowStyle = 'modal';
           STR={'\fontsize{11}Custom Honeycomb(HC) lattice cross-sections: ';
                'This UI provides five different customized cross-sections for HC cross-sections.';
                'Please use this UI \bfwith the Sketch tab to convert lines into bundles\rm as Top-down approach. After conversion, these bundles can be saved into the library for Bottom-up approach.';
                '';
                '(1): Select one of five custom cross-sections to change with the popup menu.';
                '(2): Use left clicks on the centers of the circles(cylinders)(\color{blue}blue dot\color{black}) to make it into \color{red}red\color{black} as selected.';
                '(3): Right clicks can do reversely.';
                '(4): The bottom-left text shows the numbers of two groups of cylinders which have opposite directions.';
                '  These two numbers have to the same due to the scaffold algorithm. '
                '(5): Once finished, press \bf''Enter''\rm to update the data for converting lines into bundles with the table in the Sketch tab. ';
                'For cross-sections which can''t be paired, the program will send a warning and the data won''t be updated. ';
                '';
                
                '\bf[c]\rm: clear all selected cylinders.';
                
                '';
                ' (The font size of text can be assigned by clicking on the text. The font size of the popup menu will follow after clicking it.)';
             
           
              
               
               } ;
                                 
                      f = msgbox(STR ,'Instructions', 'help' ,CreateStruct);    

       end

% \color[rgb]{0,0.5,0}


end

