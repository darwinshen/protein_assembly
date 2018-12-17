%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%     PRINT PREDICATES: print the output in pdb format  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% stampa/2
%%% INPUT: a Stream or the string 'video' and a String
%%% EFFECT: It prints the String in the stream (a file)
%%%         or in the video
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 13/11/2009

stampa(video,String) :-
    !, write(String).
stampa(Stream,String) :-
    write(Stream,String).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% print_results/4 save the results of a computation
%%% (namely, the Primary, and Tertiary structures of the protein
%%% together with some statistics, like running time and
%%% the list Src of the protein from where fragmnets are
%%% selected) on a "pdb" file.
%%% In the "pdb" format spaces matter.
%%%%%%%%%%
%%% auxiliary: print_results/5, print_results/6 that initialize
%%% the print stage and
%%% print_results_rec that enters in the Primary/Tertiary/Src lists
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 11/02/2011

print_results(Primary,  Tertiary, Time, Src, Energies, Codes) :-
    print_results('protein', Primary, Tertiary, Time, Src, Energies,Codes).

print_results(ID, Primary, Tertiary, Time, Src, Energies, Codes) :-
    name(ID,ProtL),
    ct(N), name(N,Count),
    Energies=[EN|_],
    (EN < 0 -> EEN is -EN; true -> EEN=0), name(EEN,EName),
    append(["results/Output_",ProtL,"_",Count,"_",EName,".pdb"],NameCode),
    name(FilePath,NameCode), %%% Opt: open(FilePath,append,File),
    open(FilePath,write,File),
    print_results1(File, ID, Time, Primary, Tertiary,Src,Energies),
    stampa(File,primary(Primary)),stampa(File,'.\n'),
    Tertiary=[CA,_],
    stampa(File,tertiary(CA)),stampa(File,'.\n'),
    stampa(File,codes(Codes)),stampa(File,'.\n'),
    close(File),!. %%%  opt: add print_results(video,ID,Time,Primary,Tertiary,Src).

print_results1(MEDIA, ID, Time, Primary, Tertiary,Src,Energies) :-
     Energies = [_,TorsionEnergy,ContactEnergy,TorsContributions,CorrContributions,ORContrib],   
     stampa(MEDIA,'HEADER      '), stampa(MEDIA,ID), stampa(MEDIA,'\n'),
     stampa(MEDIA,'REMARK      '), stampa(MEDIA,'Solution computed in time: '), stampa(MEDIA,Time), stampa(MEDIA,' s.   '), stampa(MEDIA,'\n'),
     write_remarks(1,MEDIA,TorsionEnergy,ContactEnergy,TorsContributions,CorrContributions,ORContrib),
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     stampa(MEDIA,'MODEL       1'),
     stampa(MEDIA,'\n'),
     print_results_rec(MEDIA, 1 , Primary, Tertiary,Src ),
     print_conect(MEDIA,1,Primary),
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     stampa(MEDIA,'ENDMDL'),stampa(MEDIA,'\n'),
     stampa(MEDIA,'END'),   stampa(MEDIA,'\n').

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% print_all_models/2 prints all the computed models,
%%%    stored in List, in a unique pdb file
%%% aux: print_all_models/3 and print_all/3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Revised 11/02/2011

print_all_models(ID):-
    findall(X,s(X),List),
    name(ID,ProtL),
    append(["results/Output_",ProtL,"_all.pdb"],NameCode),
    name(FilePath,NameCode),
    open(FilePath,write,File), 
    print_all_models(File, ID, List),
    close(File), !.

print_all_models(MEDIA, ID, List) :-
     stampa(MEDIA,'HEADER      '), stampa(MEDIA,ID),  stampa(MEDIA,'\n'),
     write_remarks(0,MEDIA,0,0,0,0,0),
     print_all(MEDIA,List,1),
     stampa(MEDIA,'END'), stampa(MEDIA,'\n').

print_all(_MEDIA,[],_N).
print_all(MEDIA,[S|Sol],N):-
     (S=[Primary,Tertiary,Src,PCost,Time,[Energy,TorsionEnergy,ContactEnergy,TorsContributions,CorrContributions,OR]],!;
%%% AGO: non capito quando si puo' presentare questa alternativa.
%%% per ora l'ho lasciata.
      Primary=[a,a,a,a],  Tertiary=S, Src=[_], PCost=0, Energy=0),

     stampa(MEDIA,'REMARK      '),stampa(MEDIA,'Solution computed in time: '),stampa(MEDIA,Time),stampa(MEDIA,' s.   '),stampa(MEDIA,'\n'),
     stampa(MEDIA,'REMARK      '),stampa(MEDIA,'Solution (log of) Probability: '),stampa(MEDIA,PCost),stampa(MEDIA,'\n'),
     stampa(MEDIA,'REMARK      '),stampa(MEDIA,'Solution Energy: '),stampa(MEDIA,Energy),stampa(MEDIA,'\n'),
     stampa(MEDIA,'REMARK      '),stampa(MEDIA,'ETors: '),stampa(MEDIA,TorsionEnergy),stampa(MEDIA,'\n'),
     stampa(MEDIA,'REMARK      '),stampa(MEDIA,'EContrTors: '),stampa(MEDIA,TorsContributions),stampa(MEDIA,'\n'),
     stampa(MEDIA,'REMARK      '),stampa(MEDIA,'EContrCorr: '),stampa(MEDIA,CorrContributions),stampa(MEDIA,'\n'),
     stampa(MEDIA,'REMARK      '),stampa(MEDIA,'EOrientation: '),stampa(MEDIA,OR),stampa(MEDIA,'\n'),
     stampa(MEDIA,'REMARK      '),stampa(MEDIA,'ECont: '),stampa(MEDIA,ContactEnergy),stampa(MEDIA,'\n'),
     stampa(MEDIA,'MODEL      '), stampa(MEDIA,N), stampa(MEDIA,'\n'),
     print_results_rec(MEDIA, 1 , Primary, Tertiary,Src ),
     stampa(MEDIA,'ENDMDL'), stampa(MEDIA,'\n'),
     N1 is N+1,
     print_all(MEDIA,Sol,N1).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% print_results_rec prints the lines corresponding to
%%% the positions of the CA atoms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%17/11/2009

print_conect(_MEDIA,_I,[_]).
print_conect(MEDIA,I,[_,A|R]):-
    (
     options(use_centroids),
     I1 is 2*I,
     I2 is 2*I+1,
    stampa(MEDIA,'CONECT  '), %%% 1--7
    (I1>99;I1>9,stampa(MEDIA,' ');I1>0,stampa(MEDIA,'  ');true),
    stampa(MEDIA,I1), %% 7-8 (atom number)
    stampa(MEDIA,'  '),
    (I2>99;I2>9,stampa(MEDIA,' ');I2>0,stampa(MEDIA,'  ');true),
    stampa(MEDIA,I2), %% 7-8 (atom number)
    stampa(MEDIA,'\n'),
     I3 is 2*I+2;
     \+options(use_centroids),
     I1 is I,
     I3 is I+1
    ),

    stampa(MEDIA,'CONECT  '), %%% 1--7
    (I1>99;I1>9,stampa(MEDIA,' ');I1>0,stampa(MEDIA,'  ');true),
    stampa(MEDIA,I1), %% 7-8 (atom number)
    stampa(MEDIA,'  '),
    (I3>99;I3>9,stampa(MEDIA,' ');I3>0,stampa(MEDIA,'  ');true),
    stampa(MEDIA,I3), %% 7-8 (atom number)
    stampa(MEDIA,'\n'),
    IN is I+1,
  print_conect(MEDIA,IN,[A|R]).

print_results_rec(MEDIA, I , [Amino|R], [[X,Y,Z|N],[XC,YC,ZC|NC]],SRC) :-

    (
     options(use_centroids),
     I1 is 2*I,
     I2 is 2*I+1;
     \+options(use_centroids),
     I1 is I,
     I2 is I
    ),

    reference_aa_input_prot(Ofs),
    AANumber is I-1+Ofs,
    !,
    (SRC = []        -> ID = 'NULL', SRC1 = [];
     SRC = [ID|SRC1] -> true),
    stampa(MEDIA,'ATOM    '), %%% 1--6
    (I1>99;I1>9,stampa(MEDIA,' ');I1>0,stampa(MEDIA,'  ');true),
    stampa(MEDIA,I1), %% 7-8 (atom number)
    stampa(MEDIA,'  CA  '), %%% 9-17
    name_conversion(AminoName,Amino,_), %% 18-20
    name(AminoNameN,AminoName),
    stampa(MEDIA,AminoNameN),
    stampa(MEDIA,' A '), %% 21-23
    (AANumber>99;AANumber>9,stampa(MEDIA,' ');AANumber>0,stampa(MEDIA,'  ');true),
    stampa(MEDIA,AANumber), %% 24-25  (residue number)
    stampa(MEDIA,'    '), %% 26-30
    %%%
    stampa_coordinata(MEDIA,X),
    stampa_coordinata(MEDIA,Y),
    stampa_coordinata(MEDIA,Z),
    %%%% annotation for original template source
    stampa(MEDIA,'  1.00  1.00              '),stampa(MEDIA,ID),stampa(MEDIA,'\n'),
    %%%
    (
    options(use_centroids),
    (ground([XC,YC,ZC]),
    (SRC = []        -> ID = 'NULL', SRC1 = [];
     SRC = [ID|SRC1] -> true),
    stampa(MEDIA,'ATOM    '), %%% 1--6
    (I2>99;I2>9,stampa(MEDIA,' ');I2>0,stampa(MEDIA,'  ');true),
    stampa(MEDIA,I2), %% 7-8 (atom number)
    stampa(MEDIA,'  CB  '), %%% 9-17
    name_conversion(AminoName,Amino,_), %% 18-20
    name(AminoNameN,AminoName),
    stampa(MEDIA,AminoNameN),
    stampa(MEDIA,' A '), %% 21-23
    (AANumber>99;AANumber>9,stampa(MEDIA,' ');AANumber>0,stampa(MEDIA,'  ');true),
    stampa(MEDIA,AANumber), %% 24-25  (residue number)
    stampa(MEDIA,'    '), %% 26-30
    %%%
    stampa_coordinata(MEDIA,XC),
    stampa_coordinata(MEDIA,YC),
    stampa_coordinata(MEDIA,ZC),
    %%%% annotation for original template source
    stampa(MEDIA,'  1.00  1.00              '),stampa(MEDIA,ID),stampa(MEDIA,'\n')
    %%%
    ;
    \+options(use_centroids),
    true)
    ;
    true
    ),

    IN is I + 1,
     print_results_rec(MEDIA, IN, R, [N,NC],SRC1).
print_results_rec(_, _ , [], [[],[]],_).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% stampa_coordinata/2 prints in the correct scale and with the
%%% proper identation the passed coordinate S in the pdb
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 18/11/2009

stampa_coordinata(MEDIA,S) :-
     resolution(Res),
     T is integer(S*1000/Res),
     U is T/1000,
     add_zeros_before(MEDIA,T,U),
     stampa(MEDIA,U),
     add_zeros_after(MEDIA,T).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% add_zero_before allows right alignment of successive symbol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 18/11/2009

add_zeros_before(MEDIA,X,Y) :-
     ( Y>99.999,!, stampa(MEDIA,' ');
       Y>9.999,!, stampa(MEDIA,'  ');
       X>=0,!, stampa(MEDIA,'   ');
       Y< -99.999,!;
       Y< -9.999,!,stampa(MEDIA,' ');
       X<0,!,stampa(MEDIA,'  ');
       true). %%% ?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% add_zero_after allows alignment of successive symbol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 18/11/2009

add_zeros_after(MEDIA,X) :-
     0 is X mod 1000, !, stampa(MEDIA,'00').
add_zeros_after(MEDIA,X) :-
     0 is X mod 100,  !, stampa(MEDIA,'00').
add_zeros_after(MEDIA,X) :-
     0 is X mod 10,  !,  stampa(MEDIA,'0').
add_zeros_after(_MEDIA,X) :-
     X mod 10 > 0.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% write_remarks prints information about the protein and
%%%      search parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 17/11/2009

write_remarks(IND,MEDIA,TorsionEnergy,ContactEnergy,TorsContributions,CorrContributions,ORContrib) :-
     min_distance(Min_d),
     min_distance_cg(Min_dcg),
     max_distance(Max_d),
     distanceSsbond(DistSSbond),
     distanceCaCa(DistCACA),
     distanceCaCaStep2(DistCACAs2),
     next_CaCadistance(NextCACA),
     options(max_most_frequent(MMF)),
     options(max_non_optimal_choices(MNOC)),
     resolution(Res),
     %%%% Comment lines
     stampa(MEDIA,'REMARK      '),stampa(MEDIA,'min_distance('),stampa(MEDIA,Min_d),stampa(MEDIA,').\n'),
     stampa(MEDIA,'REMARK      '),stampa(MEDIA,'min_distance_cg('),stampa(MEDIA,Min_dcg),stampa(MEDIA,').\n'),
     stampa(MEDIA,'REMARK      '),stampa(MEDIA,'max_distance('),stampa(MEDIA,Max_d),stampa(MEDIA,').\n'),
     stampa(MEDIA,'REMARK      '),stampa(MEDIA,'distanceSsbond('),stampa(MEDIA,DistSSbond),stampa(MEDIA,').\n'),
     stampa(MEDIA,'REMARK      '),stampa(MEDIA,'distanceCaCa('),stampa(MEDIA,DistCACA),stampa(MEDIA,').\n'),
     stampa(MEDIA,'REMARK      '),stampa(MEDIA,'distanceCaCaStep2('),stampa(MEDIA,DistCACAs2),stampa(MEDIA,').\n'),
     stampa(MEDIA,'REMARK      '),stampa(MEDIA,'next_CaCadistance('),stampa(MEDIA,NextCACA),stampa(MEDIA,').\n'),
     stampa(MEDIA,'REMARK      '),stampa(MEDIA,'resolution('),stampa(MEDIA,Res),stampa(MEDIA,').\n'),
     stampa(MEDIA,'REMARK      '),stampa(MEDIA,'max_most_frequent('),stampa(MEDIA,MMF),stampa(MEDIA,').\n'),
     stampa(MEDIA,'REMARK      '),stampa(MEDIA,'max_non_optimal_choices('),stampa(MEDIA,MNOC),stampa(MEDIA,').\n'),
     stampa(MEDIA,'REMARK      \n'),
     best(_BestEnergy),
%%% ATTENTO CON BEST ENERGY nel caso LS. 
     (IND = 0 -> true;
      IND > 0 ->
      (NewEnergy is ContactEnergy + TorsContributions + CorrContributions + ORContrib,    
       stampa(MEDIA,'REMARK      '),stampa(MEDIA,'Total Energy: '),stampa(MEDIA,NewEnergy),stampa(MEDIA,'\n'),
       stampa(MEDIA,'REMARK      '),stampa(MEDIA,'          ETors: '),stampa(MEDIA,TorsionEnergy),stampa(MEDIA,'\n'),
       stampa(MEDIA,'REMARK      '),stampa(MEDIA,'     EContrTors: '),stampa(MEDIA,TorsContributions),stampa(MEDIA,'\n'),
       stampa(MEDIA,'REMARK      '),stampa(MEDIA,'     EContrCorr: '),stampa(MEDIA,CorrContributions),stampa(MEDIA,'\n'),
       stampa(MEDIA,'REMARK      '),stampa(MEDIA,'          ECont: '),stampa(MEDIA,ContactEnergy),stampa(MEDIA,'\n'),
       stampa(MEDIA,'REMARK      '),stampa(MEDIA,'   Eorientation: '),stampa(MEDIA,ORContrib),stampa(MEDIA,'\n'))).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% nl_err reports an error during loading
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 17/11/2009

nl_err:- write(user_error,'\n').

%% computes the two inner centroids
%% for each aminoacid in the class

writeTemplateWithCentroids(Term):-
      Term=tupla([C1,C2,C3,C4],[X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4],Prob,Id,ProtName),
      get_list_centroids(C2,[X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3],ListCentroids2),
      get_list_centroids(C3,[X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4],ListCentroids3),
      write('tuple('),
      write([C1,C2,C3,C4]),write(','),
      write([X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4]),write(','),
      write([ListCentroids2,ListCentroids3]),write(','),
      write(Prob),write(','),
      write(Id),write(','''),
      write(ProtName),
      write(''').\n').

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

write_dom_ter([]).   
write_dom_ter([X,Y,Z|R]):-
  fd_set(X,SetX),  
  fd_set(Y,SetY),  
  fd_set(Z,SetZ),  
  write(p(SetX,SetY,SetZ)),nl,
  %write(p(ListX,ListY,ListZ)),nl,
  write_dom_ter(R).

write_dom_list_mat([]). %% special for distances
write_dom_list_mat([[N,L]|R]):-
  writepadding(N),
  write_dom_list(L),nl,
  write_dom_list_mat(R).

write_dom_list_mat_en([]).  %% special for energy
write_dom_list_mat_en([[N,L]|R]):-
  writepadding(N),
  write_dom_list_en(L),nl,
  write_dom_list_mat_en(R).

writepadding(0).
writepadding(N):-
  N>0,!,
  write('   '),
  N1 is N-1,
  writepadding(N1).
  

write_dom_list([]).   
write_dom_list([X|R]):- 
  fd_max(X,SetX),    %% mostra massimo  
  fd_size(X,SetSizeX),  
  resolution(Res),   
  (SetX=sup,!,
   write('oo');
   %Scale is integer(sqrt(SetX )/ Res),
   Scale is sqrt(SetX )/ Res,
   (Scale<10,!,write(' ');true),
   write(Scale)
   ),  
   %% separator gives idea of the domain size
   (SetSizeX=sup,!,write('U');
    SetSize is integer(sqrt(SetSizeX )/ Res),
    ( SetSize =< 1,!,
     write('.');
     SetSize=< 10,!,
     write('o');
     SetSize =< 100,!,
     write('O');     
     write('*')
    )
    ),
%  fd_min(X,Min),
%  fd_max(X,Max),
%  (Max \= sup,!,
%  Delta is Max-Min,
%  write(Delta),write(' ');
%  write('Inf ')),
  write_dom_list(R).



write_dom_list_en([]).   
write_dom_list_en([X|R]):-
  fd_set(X,Set),  
  write(Set),
  %fd_max(X,SetX),  
  fd_size(X,SetSizeX),  
%  resolution(Res),   
%  (SetX=sup,!,
%   write('oo');
%   Scale is integer(SetX / 1000),
%   (Scale<10,!,write(' ');true),
%   write(Scale)
%   ),  
   %% separator gives idea of the domain size
   (SetSizeX=sup,!,write('U');
    SetSize is integer(SetSizeX/1000),
    ( SetSize =< 1,!,
     write('.');
     SetSize=< 10,!,
     write('o');
     SetSize =< 100,!,
     write('O');     
     write('*')
    )
    ),
%  fd_min(X,Min),
%  fd_max(X,Max),
%  (Max \= sup,!,
%  Delta is Max-Min,
%  write(Delta),write(' ');
%  write('Inf ')),
  write_dom_list_en(R).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%     END        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

