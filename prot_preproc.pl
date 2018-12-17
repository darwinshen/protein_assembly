%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% prot preproc retrieves the protein data from the DB
%%%% and generates the desired input data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%       DATABASE HANDLING      %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Data base handling: load only useful templates
%%%                     and create a local db to
%%%                     be compiled (more efficient)
%%% the db includes the list of templates for each 4-ple
%%% and the compatible next concatenation together
%%% with the rotation matrix
%%%
%%% Assuming that the template data is loaded,
%%% generates the energy contributions
%%% for torsional angles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% reset_db forces the recomputation of the templates db
%%% and of the associated energy contributions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%12/11/2009

%%% The DB for the specific protein is added with the centroids
%%% for positions 2 and 3 (over ca1,ca2,ca3,ca4) as a pair
%%% of lists of centroids positions for each specific amino acid

reset_db(ID):-
   %%build filename
    name(ID,LID),
    append(["temp/",LID,"-current.db.pl"],Name),
    name(FileName,Name),
    on_exception(_, (
                   delete_file(FileName),
                   write('Database successfully reset for protein '),write(ID)
                  ),
   write('The Database was already reset')).

init_data(ID):-
   load_db(ID),
   load_torsional_energy(ID).

load_db(ID):-
   proteina(ID,Primary,Secondary,_),
   seen,
   told,
   %%build filename
    name(ID,LID),
    append(["temp/",LID,"-current.db.pl"],Name),
    name(FileName,Name),
   %%% is the current db correct (same protein sequence described?)
   (
     on_exception(_,see(FileName),
     (write('The database for the protein not cached\n'),fail)),
     read(sequence(Primary)),!,
     seen    %%% can compile
     ;
     write('creating template db for protein '),write(ID),nl,
     seen,
     see('db/templates.pl'),
     tell(FileName),
     write(sequence(Primary)),write('.\n'),
     retractall(temp(_)), %% temporarily asserted during creation
     read(A),
%    trace,
     loop(A,Primary,Secondary),
     seen,
     told,
     retractall(temp(_)), %% temporarily asserted during creation
     write('current db created'),nl
   ),
   write('compiling template db for protein '),write(ID),nl,
   abolish(sequence/1,[force(true)]),
   abolish(class/2,[force(true)]),
   abolish(tuple/6,[force(true)]),
   abolish(next/3,[force(true)]),
   abolish(originalCodes/2,[force(true)]),

   compile(FileName),
   write('compile ok'),nl.

loop(Term,Primary,Secondary):-
   (Term=class(_,_),!,
    %write(user_error,Term),nl_err,
    write(Term),write('.\n'),
    assert(temp(Term)),
    read(A),
    loop(A,Primary,Secondary);

    convert_list_temp(Primary,ConvertedPrimary),  %% uses temp facts
    write(user_error,ConvertedPrimary),nl_err,
    write(user_error,'scan templates: '),
    loop1(Term,[],[],ConvertedPrimary,Secondary)).

loop1(Term,Templ,SpecialTempl,Seq,Secondary):-
    (
     Term=tupla(List,_Points,Freq,Id,_Prot),!,

     %% assert tuple only if it's in the Seq(uence)
     (sublist(Seq,List,_,_,_),Freq>0,!,
%write(user_error,test(ok,Seq,List)),nl_err,
      %write(user_error,Term),write(user_error,'\n'),
      %write(Term),write('.\n'),
      writeTemplateWithCentroids(Term),
      assert(temp(Term)),
      NewTempl=[Id|Templ],
      NewSpecialTempl=SpecialTempl,
      write(user_error,'*')
      ;
      List=[C,C,C,C],C<0,Freq>0,!,  %% special
%write(user_error,test(sp,Seq,List)),nl_err,
      NewTempl=Templ,
      writeTemplateWithCentroids(Term),
      %write(Term),write('.\n'),
      NewSpecialTempl=[Id|SpecialTempl],
      write(user_error,'*')
      ;
%write(user_error,test(no,Seq,List)),nl_err,
      NewTempl=Templ,
      NewSpecialTempl=SpecialTempl
     ),

     read(A),
     loop1(A,NewTempl,NewSpecialTempl,Seq,Secondary); %% other loop1 test

     %%% verify whether special template needs to be introduced
     %write(user_error,SpecialTempl),nl_err,
     nl_err,
     (check_no_need_special(Seq),check_no_need_secondary(Secondary),!,
      write(user_error,'Templates contained in db'),nl_err,
      AllTempl=Templ;
      write(user_error,'Adding special general templates and SSE template'),nl_err,
      %%% if not all templates are supported -> add special general templates (assert them and add codes)
      append(Templ,SpecialTempl,AllTempl)
     ),


     Term=originalCodes(ListCodes,PName),
     write('originalCodes('), write(ListCodes), write(',\''), write(PName), write('\').\n'),
     write(user_error,fixedcodes(PName)),nl_err,
     write(user_error,codes(AllTempl)),nl_err,
     write(user_error,'Now load next: '),nl_err,
     seen,
     see('db/nextrange.pl'),
     read(range(RangeList)),
     seen,
     sort(AllTempl,AllTemplS),
     write(user_error,'Locating files with next...'),nl_err,

     setupnextfilename(AllTemplS,RangeList,FileNames,PartitionedTempl,_TG),   %% for each AllTempl code, find the corresponding file and cluster AllTempl into same group
%write(FileNames),nl,
%write(PartitionedTempl),nl,
%write(TG),nl,
     %write(FileNames),nl,
     loadnext(PartitionedTempl,FileNames,AllTemplS),
     write(user_error,ok),nl_err
     ). %% go on with next (for each code in AllTempl, load the next predicate)

get_list_centroids(Class,Points,ListCentroids):-
      (Class>=0,!,
       findall(AA,(temp(class(AA3,Class)),
                   name_conversion(AA3, AA, _)),
               ListAA)
       ;
       findall(AA,(temp(class(AA3,_)),
                   name_conversion(AA3, AA, _)),
               ListAA)
      ),
      %%for each AA in ListAA, compute the corresponding centroid
      calc_centroid_list(ListAA,ListCentroids,Points).

calc_centroid_list([],[],_).
calc_centroid_list([AA|R],[[AA,Centroid]|R1],Points):-
   calc_centroid(AA,Centroid,Points),
   calc_centroid_list(R,R1,Points).

calc_centroid(g,[X2,Y2,Z2],[_X1,_Y1,_Z1,X2,Y2,Z2,_X3,_Y3,_Z3]):-!.
calc_centroid(AA,[CxND,CyND,CzND],[X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3]):-
        V1x is X2-X1, V1y is Y2-Y1, V1z is Z2-Z1,
        V2x is X3-X2, V2y is Y3-Y2, V2z is Z3-Z2,
        N1 is sqrt(V1x*V1x+V1y*V1y+V1z*V1z),
        N2 is sqrt(V2x*V2x+V2y*V2y+V2z*V2z),
        V1xN is V1x / N1, V1yN is V1y / N1, V1zN is V1z / N1,
        V2xN is V2x / N2, V2yN is V2y / N2, V2zN is V2z / N2,

        %centroid_ang1(AA,Ang1),
        centroid_ang2(AA,Ang2),
        centroid_tors(AA,Tors),

/* Sottraggo da v1 la componente lungo v2 in modo da avere v1 e v2 ortogonali */
        TempA is V1xN*V2xN+V1yN*V2yN+V1zN*V2zN,
        V1xT is V1xN - TempA*V2xN,
        V1yT is V1yN - TempA*V2yN,
        V1zT is V1zN - TempA*V2zN,
        NT is sqrt(V1xT*V1xT+V1yT*V1yT+V1zT*V1zT),
        V1xTN is V1xT / NT, V1yTN is V1yT / NT, V1zTN is V1zT / NT,

        prod_vect(V1xTN,V1yTN,V1zTN,V2xN,V2yN,V2zN,V3xN,V3yN,V3zN),

        B0 is cos(Ang2/180*3.1415),
        B1 is sin(Tors/180*3.1415)* sqrt(1-B0*B0),
        B2 is cos(Tors/180*3.1415)* sqrt(1-B0*B0),

        R00 is V2xN, R01 is V2yN, R02 is V2zN,
        R10 is V3xN, R11 is V3yN, R12 is V3zN,
        R20 is -V1xTN, R21 is -V1yTN, R22 is -V1zTN,

        D is
        R00 * R11 * R22 +
        R01 * R12 * R20 +
        R02 * R10 * R21 -
        R02 * R11 * R20 -
        R01 * R10 * R22 -
        R00 * R12 * R21,

        Dx is
          B0 * (R11 * R22 - R12 * R21)
        + B1 * (R21 * R02 - R22 * R01)
        + B2 * (R01 * R12 - R02 * R11),

        Dy is
          B0 * (R12 * R20 - R10 * R22)
        + B1 * (R22 * R00 - R20 * R02)
        + B2 * (R02 * R10 - R00 * R12),

        Dz is
          B0 * (R10 * R21 - R11 * R20)
        + B1 * (R20 * R01 - R21 * R00)
        + B2 * (R00 * R11 - R01 * R10),

        CxV is Dx/D, CyV is Dy/D, CzV is Dz/D,

        centroid_dist(AA,Dist),
        %% normalize
        %Norm is sqrt(CxV*CxV+CyV*CyV+CzV*CzV),
        CxND is integer(1000*(CxV*Dist+X2))/1000.0,
        CyND is integer(1000*(CyV*Dist+Y2))/1000.0,
        CzND is integer(1000*(CzV*Dist+Z2))/1000.0.



%%% true if all 4-tuples in the sequence are covered by asserted facts
check_no_need_special([]).
check_no_need_special([_]).
check_no_need_special([_,_]).
check_no_need_special([_,_,_]).
check_no_need_special([A,B,C,D|R]):-
%write(user_error,test([A,B,C,D])),nl,
   temp(tupla([A,B,C,D],_,_,_,_)),!,
   check_no_need_special([B,C,D|R]).

check_no_need_secondary(Secondary):-
     %%% if secondary structure is present -> force the presence for the best generic alpha and beta pattern
     nonmember(helix(_,_),Secondary),
     nonmember(strand(_,_),Secondary),!.

%%% for each template open the correct file and read the next associated
%%% Group is the last range processed. if C in the the same range,
%%% no newFilePath is added and C is only inserted in Part
setupnextfilename([],_,[],[],[]).
setupnextfilename([C|Templ],RangeList,NewFilePath,NewPart,NewGroup):-
  %write(user_error,C),nl_err,
  setupnextfilename(Templ,RangeList,FilePath,Part,Group),
  %write(user_error,s(C,FilePath,Part,Group)),nl_err,

  (
    Group=[[A-B]|_],
    A=<C,C=<B,!,    %%% can recycle last one
    NewFilePath=FilePath,
    Part=[CurrentPart|R],
    NewPart=[[C|CurrentPart]|R],
    NewGroup=Group
    ;
       %%%% create new
    findall(NA-NB,(member([NA,NB],RangeList),NA=<C,C=<NB),[RA-RB]),
    name(RA,AL),
    name(RB,BL),
    append("db/next-",AL,S1),
    append(S1,"-",S2),
    append(S2,BL,S3),
    append(S3,".pl",S4),
    name(FP,S4),
    NewFilePath=[FP|FilePath],
    NewPart=[[C]|Part],        %% start new partition
    NewGroup=[[RA-RB]|Group]
  ).

loadnext([],[],_).
loadnext([SetT|Templ],[FN|R],TemplAll):-
  see(FN),
  read(Data),
  write(user_error,loadFrom(FN,SetT)),
  loop2(SetT,Data,TemplAll),
  nl_err,
  seen,
  loadnext(Templ,R,TemplAll).

loop2(SetT,Term,Templ):-     %% Templ is the list of codes that are extracted from the db based on the protein sequence
                             %% SetT contains the available templates in current file
    Term=next(_,_,_),!,    %% not end_of_file
    (Term=next(Id1,Id2,_Matrix),
     member(Id1,SetT),
     member(Id2,Templ),!,
     %write(nextOK(SetT,Id2)),nl,
     write(user_error,'*'),
     write(Term),write('.\n');
     true),
    read(A),
    loop2(SetT,A,Templ).
loop2(_,_Term,_Templ). %% last failed (ok)

%%%%%% creates a pl file from the database, with specific torsinal data for the protein
%%%%%% The content is:
%%%%%% tors(C1,C2,List).     % for each class C1 and C2, the list of pairs [Torsional bin,Energy]
%%%%%% corr(List).           % the list of correlation triples [T1,T2,Energy] for two specific consecutive torsions
%%%%%% tors_templ(ID,Tors).  % the specific torsion for each template (ready to be used with previous facts)

load_torsional_energy(ID):-
   seen,
   told,
   %%build filename
    name(ID,LID),
    append(["temp/",LID,"-currentTorsdb.pl"],Name),
    name(FileName,Name),
   sequence(Primary),
   %%% is the current torsional db correct (same protein sequence described?)
   (
     on_exception(_,see(FileName),(write('The database for torsional energy not cached\n'),fail)),
     read(sequenceT(Primary)),!,
     seen    %%% can compile
     ;
     write('creating current torsional energy db'),nl,
     seen,
     load_torsional_data,

     extract_torsional_angles,

     findall(C,temp(C),LAll),

     write('Writing torsions to file...\n'),
     tell(FileName),

     write(sequenceT(Primary)),write('.\n'),
     write_list_as_prolog_facts(LAll),

     told,
     write('Writing torsions to file ok\n'),
     retractall(temp(_)), %% temporarily asserted during creation
     write('current torsional db created'),nl
   ),
   write('compiling current torsional energy db'),nl,
   abolish(sequenceT/1,[force(true)]),
   abolish(prof/3,[force(true)]),
   abolish(corr/1,[force(true)]),
   abolish(tors_templ/2,[force(true)]),
   compile(FileName),
   write('compile ok'),nl.


extract_torsional_angles:-
   findall(ID,tuple(_,_,_,_,ID,_),List), %% from current db get list of template codes
   %write(List),nl,
   associate_torsional_energy(List)
   %,write(ok),nl
   .

associate_torsional_energy([]).
associate_torsional_energy([ID|R]):-
  tuple(_,Templ,_Centroids,_,ID,_),
  get_tors_angle(Templ,Ang),
  assert(temp(tors_templ(ID,Ang))),
  associate_torsional_energy(R).

get_tors_angle([X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4],Ang1):-
  V1x is X2-X1,  V1y is Y2-Y1,  V1z is Z2-Z1,
  V2x is X3-X2,  V2y is Y3-Y2,  V2z is Z3-Z2,
  V3x is X4-X3,  V3y is Y4-Y3,  V3z is Z4-Z3,

  prod_vect(V2x,V2y,V2z,V1x,V1y,V1z,Ix,Iy,Iz),
  prod_vect(V2x,V2y,V2z,Ix,Iy,Iz,Jx,Jy,Jz),
  normalize(Ix,Iy,Iz,IxN,IyN,IzN),
  normalize(Jx,Jy,Jz,JxN,JyN,JzN),
  normalize(V3x,V3y,V3z,V3xN,V3yN,V3zN),
  prod_scal(IxN,IyN,IzN,V3xN,V3yN,V3zN,PrI),
  prod_scal(JxN,JyN,JzN,V3xN,V3yN,V3zN,PrJ),

  Ang is atan2(-PrI,PrJ)/atan2(1,0)/2*180,

  Ang1 is (integer(floor((Ang+3)/5)))*5-3     %% compute bin
  .

%% assert profile and correlation as temp(prof(C1,C2,Data)) and temp(corr(List))
%%   Data is a list of pairs (Torsional angle, energy)
%%   List is a list of triples (T angle1, T angle2, energy)
%%   One class for each amino acid

load_torsional_data:- 
  write('Generation of torsional energy contributions.'),nl,
  %findall(C,class(_,C),L),
  %remove_dups(L,L1),
  
  findall(X,(name_conversion(X1,_,_),name(X,X1)),L1),
  
  length(L1,Len),
  %write(L1),
  write('There are '),write(Len),write(' aa classes.'),nl,
  retractall(temp(_)),
  write('Loading profiles...\n'),
  rec_load_profile(L1,L1),
  write('Loading profiles done\n'),
  read_torsion_correlation.

rec_load_profile([],_L1).
rec_load_profile([C|R],L1):-
  rec_load_profile1(C,L1),
  rec_load_profile(R,L1).

rec_load_profile1(_C1,[]).
rec_load_profile1(C1,[C2|R]):-
  read_torsion_profile_class(C1,C2,L),
  
  %% convert from 3 letters to 1 letter
  
  name(C1,C1s),
  name_conversion(C1s,C1l,_),
  name(C2,C2s),
  name_conversion(C2s,C2l,_),
  
  assert(temp(prof(C1l,C2l,L))),
  rec_load_profile1(C1,R).

%%% torsional energy handling

%%% returns the list of torsional contribution for pair of central classes C1 and C2 (position 1,2 in the template 0,1,2,3)
%%% the energy contributions are multiplied by 1000 for FD
read_torsion_profile_class(C1,C2,L):-
  seen,

  name(C1,NC1),
  name(C2,NC2),
  append("pmf/CA1_CA2_CA3_CA4_clean_",NC1,T1),
  append(T1,"_",T2),
  append(T2,NC2,T3),
  append(T3,".pmf",FNL),
  name(FN,FNL),
  %write(FN),nl,
  see(FN),

  collect_pairs(L),
  %write(L),
  seen.

collect_pairs(Out):-
  (
  get_number(A,11),
  get_number(B,12),
  A1 is integer(A-0.5),      %%the bins are shifted by half degree 
  B1 is integer(B*1000),
  get_char(_), %%% new line
  !,collect_pairs(R),
  Out=[[A1,B1]|R]
  ;
  Out=[]
  ).

%%% L is a list of triples (ang, ang, energy) of correlation between torsional angles
read_torsion_correlation:-
  seen,
  see('pmf/t1_t2_new.pmf'),
  write('Loading correlations...'),nl,
  collect_triples(L),
  write('Loading correlations done\n'),
  %write(L),
  assert(temp(corr(L))),
  seen.

collect_triples(Out):-
  (
  get_number(A,11),
  get_number(B,12),
  get_number(C,12),
  A1 is integer(A-0.5), %%the bins are shifted by half degree 
  B1 is integer(B-0.5), %%the bins are shifted by half degree 
  C1 is integer(C*1000),
  get_char(_), %%% new line
  !,collect_triples(R),
  Out=[[A1,B1,C1]|R]
  ;
  Out=[]
  ).

get_number(A,N):-
  get_number1(A1,N),
  number_chars(A,A1).

get_number1([],0).
get_number1([C|A],N):-
  on_exception(_, get_char(C),fail),
  N1 is N-1,
  get_number1(A,N1).
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% write_list_as_prolog_facts "prints" as prolog facts
%%%  a list of atoms.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 12/11/2009

write_list_as_prolog_facts([]).
write_list_as_prolog_facts([A|L]):-
   write(A),write('.'),nl,
   write_list_as_prolog_facts(L).
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%     END        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

