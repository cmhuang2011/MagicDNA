function [OneOverHang,stapletext0]=getOverSeq( overhang_length,OtherOverhangs)

%% Josh Johnson 
% This program uses a give scaffold and a text list of other
% sequences to determine which overhang sequences will be the best. The
% basic idea is that instead of randomly generating sequences and checking
% them, this program finds the set of all possible short (usually 5 bases)
% sequences then determines how many times they appear in the scaffold.
% This creates a matrix, like an occuarance score card, that is used for
% building the rest of the overhang. It starts with the sequences that
% occur the least. It then builds the overhang by looking at the last four
% bases and using the matrix of scores to determine which new base will
% make the resulting sequence have the lowest occurance score. Since all
% sequence searching was already done to building the occurance score, this
% process is very fast. Disclaimer: this code is still work in progress

%Change to function, Chao-Min Huang
% clear all
% clc
% tic

%% Specify overhang properties here
% overhang_length = 12; % Specify your desire sequence length here  %--------------Var
gccontent = [20 80]; %specifiy minimum and maximum gc content
overhang_length=overhang_length-1;

% Import the scaffold sequence here

% fileID = fopen('p8064.txt');
% C = textscan(fileID,'%s');
% fclose(fileID);
% scaffoldtext = char(C{1});

scaffoldtext=p7560 ;


% Import your staple list here (as a .txt file)
% filename=uigetfile('*.txt','Choose scaffold file');
% SC_file=importdata(filename);

stapletext0 = importdata('ssScaf.txt','%s');  %------need to update in outer loop---Var ,cell array
stapletext0= [stapletext0 ; OtherOverhangs]    %CM

stapletext = repmat(stapletext0,5,1);
stapletext{numel(stapletext)+1}=scaffoldtext;
originalLength = numel(stapletext);
%% STEP 1: find all starting sequences that are 5 bases long
bases = categorical({'A' 'T' 'G' 'C'});

len = 7; %the length of starting sequences %------------------
all = bases;
for i=1:len-1
    temp = all.*bases';
    all = reshape(temp,[1 numel(temp)]);
end
%% This just converts the array type to something more easy to use later

seqs = cell([1 numel(all)]);
for i=1:numel(all)
     seqs{i}=strrep(char(all(i)),' ','');
end
% seqs{1}
%% Evaluate all sequences for the number of times they are found in scaffold
primercount = zeros([1, numel(all)]);
scores = cell(1,numel(all));
for s=1:numel(stapletext)
    for i=1:numel(all)
        temp = strfind(stapletext{s}, seqs{i});
        primercount(i) = primercount(i)+numel(temp);
    end
end
% clc
k =find((primercount<=prctile(primercount,10))&(primercount>-1)); %take only the sequences that are in the 10nth percentile in terms of total number of times they are found
goodprimers = seqs(k);
figure(12313);
plot(primercount)
numel(k);
%% Make a correlation matrix between sequence and number of occurances
tmatrix = containers.Map(seqs,primercount);

%% Find the best overhangs
overhangs = goodprimers;

for i=numel(goodprimers{1}):overhang_length
    overhangs = nextbestseq(overhangs,tmatrix,len); %see function notes
end

% warning('off', 'bioinfo:oligoprop:SeqLengthTooShort')
%% Check GC content
gccont = zeros([1 length(overhangs)]);
for i=1:length(overhangs)
    temp = oligoprop(char(overhangs(i)));
    gccont(i) = temp.GC;
end
k = find((gccont<gccontent(2))&(gccont>gccontent(1)));
finalists1=overhangs(k);


%% Remove strands with hairpins or dimers
hairpins = zeros([1 length(finalists1)]);
for i=1:length(finalists1)
    temp = oligoprop(char(finalists1(i)));
    hairpins(i) = size(temp.Hairpins,1);
end
h=find(hairpins<=3);
finalists2=finalists1(h);
dimers = zeros([1 length(finalists2)]);
for i=1:length(finalists2)
    temp = oligoprop(char(finalists2(i)),'Dimerlength',4);
    dimers(i) = size(temp.Dimers,1);
end
d=find(dimers<=3);
finalists3 = finalists2(d);



%% Evaluate the overhangs and return the ones with the lowest overall complementarity to the scaffold and other sequences

ohcount = zeros([1, numel(finalists3)]);
scores = cell(1,numel(finalists3));
for i=1:numel(finalists3)
    oh = finalists3{i};
    for j=1:numel(oh)-(len-1)
        ohcount(i) = ohcount(i)+tmatrix(oh(j:j+len-1));
    end
end

k = find(ohcount<=prctile(ohcount,20)) ;
bestoverhangs = finalists3(k);

% for i=1:length(bestoverhangs)
%     fprintf('\n%s', bestoverhangs{i})
% end
% toc
OneOverHang= bestoverhangs{randi(length(bestoverhangs))} ;  %random pick
end 

% %% Compare to random sequences
% 
% randos = cell(1,numel(goodprimers));
% for i=1:numel(goodprimers)
%     randos{i} = randsample('ATCG',14,true);
% end
% 
% %% Evaluate the random sequences
% 
% randocount = zeros([1, numel(randos)]);
% scores = cell(1,numel(randos));
% for i=1:numel(randos)
%     oh = randos{i};
%     for j=1:numel(oh)-(len-1)
%         randocount(i) = randocount(i)+tmatrix(oh(j:j+len-1));
%     end
% end
% 
% 
% figure
% histogram(randocount)
% hold on
% histogram(ohcount)
% 
% %% Display sequence alignments
% score = zeros(1,numel(overhangs));
% for i=1:numel(overhangs)
%     [S A] =swalign(stapletext{1},overhangs{i});
%     score(i) = S;
%     %A
%     %pause()
% end    
% 
% 
% %% Compare to random sequences
% rscore = zeros(1,numel(randos));
% for i=1:numel(randos)
%     [S A] =swalign(stapletext{1},randos{i});
%     rscore(i) = S;
%     %A
%     %pause()
% end   
% % hold off
% % figure
% % histogram(score)
% % hold on
% % histogram(rscore)
% toc

function [ seqplusone ] = nextbestseq( seq, tmatrix,len )
%nextbestseq uses a given association map (tmatrix) to determine the next best
%nucleotide to add to a set of sequences (seq). It takes the last four
%bases from each sequence and appends A,T,C, or G. Then using tmatrix it
%selects the seqeunce that is known to occur the least. If there is no
%advantage between a particular set of bases then it will pick one at
%random. It then returns the list of all sequences with the next best base
%appended to each one.

lastfour = cell(4, numel(seq));
seqplusone = cell(1, numel(seq));
nextscore = zeros(1,4);
bases = 'ATCG';
for i = 1:numel(seq)
    dummy1 = seq{i};
    for j=1:4
        lastfour{j,i} = [dummy1(end-(len-2):end) bases(j)]; %get the last 4 bases of all sequences and add one base
        nextscore(j) = tmatrix(lastfour{j,i}); %evaluate the new sequences
    end
    next = bases(find(nextscore==min(nextscore))); 
    if next>1
        next = randsample(next,1); %add a random base if two or more have then same score
    end
    seqplusone{i} = [dummy1 next]; %return new sequences
end

end