function r = getGlobalstapleLength
global stapleLength  ;


f571 = figure(571);clf ; bar(stapleLength)  ;
 f571.Name='a group of staple length' ; set(f571,'NumberTitle','off');
xlabel('staple index in this group' ) ;
ylabel('staple length ' ) ;

r = stapleLength;