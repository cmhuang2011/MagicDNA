
% script for obtaining assigned number of overhangs with custum lengths 

clear ;clc
tic
fprintf('script for obtaining assigned number of overhangs with custum lengths  \n')

% OverhangPrefer=[ 3, 8  ;5, 12 ];  % [n , lengths ]

% OverhangPrefer=[ 18, 21   ];  % [n , lengths ]
OverhangPrefer=[ 8, 25 ;10;25   ];  % [n , lengths ]

OV=[];
for k=1:size(OverhangPrefer ,1) 
    n= OverhangPrefer(k ,2) ;
    for j=1:OverhangPrefer(k ,1)
%         [k,j]
%         fprintf('state = %i %i  \n',k,j) ;
    [OneOverHang,ssDNAScaf]=getOverSeq( n, OV)      ;
%     OneOverHang
     fprintf('Seq [%i %i]  = %s  \n',k,j,OneOverHang) ;
    OV{end+1,1}=OneOverHang;
%     OV=OV ; OneOverHang ;
    end
end



fprintf('this script is over \n')
toc
return
%%
%%Get closing strand(fake scaffolds) sequences

fileIDOVandClosing = fopen('SaveSeqForFuture_BrickSeq.txt','w');

% for k=1:length(OV)
% %     SeqShort = seqrcomplement(OV{k}) ;
% %     SeqLong = seqrcomplement(OV{k+6}) ;
%     fprintf(fileIDOVandClosing,'%i %s\n',k,OV{k}) ;
% %     fprintf(fileIDOVandClosing,'%s\n',OV{k+6}) ;
% 
% %     Closing=[SeqShort SeqLong] ;
% %     duplicates= repmat(Closing,[1,6]) ;
% %     fprintf(fileIDOVandClosing,'%s\n',duplicates) ;
% %     fprintf(fileIDOVandClosing,'\n')
% end
% for k=1:6
%     SeqShort = seqrcomplement(OV{k}) ;
%     SeqLong = seqrcomplement(OV{k+6}) ;
%     fprintf(fileIDOVandClosing,'%s\n',OV{k}) ;
%     fprintf(fileIDOVandClosing,'%s\n',OV{k+6}) ;
% 
%     Closing=[SeqShort SeqLong] ;
%     duplicates= repmat(Closing,[1,6]) ;
%     fprintf(fileIDOVandClosing,'%s\n',duplicates) ;
%     fprintf(fileIDOVandClosing,'\n')
% end


for k=1:2:length(OV)
    
    fprintf(fileIDOVandClosing,'%i %s%s\n',(k+1)/2,OV{k},OV{k+1}) ;

    
end

fclose(fileIDOVandClosing);
