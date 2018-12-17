%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% By A. Dal Palu', A. Dovier, F. Fogolari, E. Pontelli
%%% July 10, 2009 --> FEBRUARY 10, 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Main Launch modes:
%%% pf(+NAME).   %%% computes the first structure (only file output) with primary structure NAME contained in the file prot-list.pl
%%% pf_enumerate(+NAME,TIME). %%% puts every structure into a different pdb file, as soon as found (timeout in seconds)
%%% pf_enumerate(+NAME). %%% puts every structure into a different pdb file, as soon as found
%%% pf_all(+NAME). %%% computes all admissible structures (and puts them as different models in the same pdb file)
%%% pf_all(+NAME,N). %%% computes the first N structures (and puts them as different models in the same pdb file)
%%% pf_lns(+NAME,TIME). %%% Performs LNS with a timeout of TIME seconds
%%%%%%%%%%
%%% add_sequence(Primary) add a new sequence addressed by "new" in the DB
%%%    We suggest however to save the structure in protlist.pl with a NEW name and use other kind of calls.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% new_options([opt1,opt2...]). %%% defines the search options
%%%   These options affect various choices in the code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

:- use_module(library(clpfd)).     :- use_module(library(terms)).
:- use_module(library(lists)).     :- use_module(library(file_systems)).
:- use_module(library(timeout)).   :- use_module(library(random)).
:- use_module(library(system)).    %:-prolog_flag(fileerrors,_,off).
:- use_module(library(aggregate)). %:-prolog_flag(unknown,_,fail).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%   MAIN PREDICATE(s) %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% pf is the main predicate.
%%% It can be called in various variants.
%%% The following 5 look for a solution per time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pf(ID) :-
   pf(ID,_Tertiary).

pf(ID,Tertiary) :-
   pf_energy(ID,Tertiary,0).

pf_energy(ID,Tertiary,InputEnergy) :-
   pf_diameter_energy(ID,Tertiary,InputEnergy,_).

pf_diameter(ID,Tertiary,Diameter) :-
   pf_diameter_energy(ID,Tertiary,_,Diameter).

pf_diameter_energy(ID,Tertiary,InputEnergy,Diameter) :-
   reset,
   (ground(InputEnergy) -> max_distance(Diameter);
    ground(Diameter)    -> InputEnergy = 0),
   pf_base(ID,Primary,Tertiary,Src,PCost,Time,Diameter,_PT,_PM, InputEnergy, Energies, Codes),
   print_results(ID,Primary,Tertiary,Time,Src, Energies,Codes),
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   format("### Total (log of) Probability = ~q~n",[PCost]),
   format("### Search Time = ~q s~n",[Time]),
   diameter(Tertiary,ED),
   format("### Effective diameter = ~3F AA~n",[ED]),
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   (options(contiguous_check) ->
       (Tertiary=[[_,_,_|CA],_CG],
        append(CA1,[_,_,_],CA),
        contiguous_check(CA1,1));
    \+options(contiguous_check) ->
        true).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% pf_enumerate  looks for all solutions, improving
%%% the Total Energy. They are saved on different pdb files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 23/11/2009

pf_enumerate(ID,Timeinseconds):-
  TOTALtime is Timeinseconds * 1000,
  time_out(pf_enumerate(ID), TOTALtime, _).

pf_enumerate(ID) :-
   reset,
   max_distance(CF),
   pf_base(ID,Primary,Tertiary,Src,_PCost,Time,CF,_PT,_PM, 0, Energies,Codes),
   diameter(Tertiary,Diam),
   retract(ct(N)), N1 is N+1, assert(ct(N1)),
   format("~a: Found structure Number ~d with Diameter ~2F\n",[ID,N1,Diam]),
   %%% PATCH DA TOGLIERE / PER SALVARE SU FILE SOLO LE MIGLIORI
   %%% Energies = [E|_], E < -18000,
   %%%%
   print_results(ID,Primary,Tertiary,Time,Src, Energies,Codes),fail.

pf_enumerate(_) :-
   ct(N), search_time(Time),
   format("### Found ~d Structures in ~1F s.\n",[N,Time]).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% pf_all looks for Howmany "improving" solutions,
%%%      saving them in the same pdb file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 17/11/2009

pf_all(ID) :-
  pf_all(ID,1000000).

pf_all(ID,Howmany) :-
   reset,
   pf_findall(ID,Howmany),
   print_all_models(ID),nl.

pf_findall(ID,Howmany):-
   max_distance(Diameter),
   pf_base(ID,Primary,Tertiary,Src,PCost,Time,Diameter,_PT,_PM,0,
        Energies,_Codes),
   assert(s([Primary,Tertiary,Src,PCost,Time,Energies])),
   ct(Ct), format("Found structure Number ~d in Time ~3F s\n",[Ct,Time]),
   (Ct<Howmany, retractall(ct(_)), Ct1 is Ct+1, assert(ct(Ct1)), fail;
    Ct>=Howmany).
pf_findall(_ID,_Howmany).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% pf_lns looks for improving solutions,
%%%      with a Timeinseconds bound (default 10s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 24/11/2009

pf_lns(ID) :-
   pf_lns(ID,10).

pf_lns(ID,Timeinseconds) :-
  pf_lns(ID,Timeinseconds,0).

pf_lns(_,_) :-
     best([0,_,_,_,_,_],_,_,_,_),!,
     write('Insufficient time for the first solution'),nl.

pf_lns(_,_) :-
   best_printing.

pf_lns(ID,Timeinseconds,Run) :-
   reset,
   TOTALtime is Timeinseconds * 1000,
   init_data(ID),
   loading_time,
   %%%%%%%%%%%%%%%%
   proteina(ID,Primary,Secondary,_),
   write(Primary),nl,
   convert_list(Primary,Codes),
   write(Codes),nl,
   %%%%%%%%%%%%%%%%%
   max_distance(Diam),
   proteina(ID,Primary,Secondary,_),
   length(Primary,N),M is N-3, good_points(1,M,Secondary,Pivots),
   %%%%%%%%%%%%%%%%%
   prepare_configuration(noOriginalTertiary,Primary,Secondary,CommandsCode,CommandsTertiary),
   originalCodes(OriginalCodes,_),


   constrain(Primary,Secondary,Tertiary,Code,PDom,Diam,_PT,_PlacedMatrices,
             Energy,TE,CE,DL,CES,1,TCP,TCC,_OriginalTertiary,_EmptyOriginalMatrices,
             OriginalCodes,CommandsCode,CommandsTertiary),
   constraint_time,
   %%%%%%%%%%%%%%%%
   retractall(best(_,_,_,_,_)),     assert(best([0,_,_,_,_,_],Primary,_,_,_)),

   retractall(last_sol(_,_,_)),     assert(last_sol(0,_,_)),
   %%%%%%%%%%%%%%%%
   name(Run,NameA),name(ID,NameB),append([NameA,"-",NameB],NameC),name(ID1,NameC),
   %%%%%%%%%%%%%%%%

   time_out(local(ID1,Code, Energy, Primary,Pivots,Tertiary, PDom, TE,CE,DL,CES,TCP,TCC),
              TOTALtime,_),
   write('*********************************************************'),nl,
   fail.


best_printing :-
     best_printing('protein').
best_printing(ID) :-
     best(Energies,Primary,Tertiary,Time,Codes),
     Energies=[Energy|_], format("### Total Energy = ~q~n", [Energy]),
     format("### Search Time = ~q s~n",[Time]),
     diameter(Tertiary,ED),
     name(ID,IDN),
     append("lns_",IDN,ID1),
     name(ID2,ID1),
     format("### Effective diameter = ~3F AA~n",[ED]),
     print_results(ID2,Primary,Tertiary,Time,[], Energies,Codes).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% pf_base is the basic predicate with all the needed parametrs
%%% used by various version of the "meta" pf predicate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pf_base(ID,Primary,Tertiary,Src,PCost,Time,CF,PlacedTemplates,PlacedMatrices,
   InputEnergy,[TOT,TorsionEnergy,ContactEnergy,
   TorsContributions,CorrContributions,OrientationEnergy],Code) :-
   init_data(ID),
   loading_time,
   write('The following options are being used: \n'),write_options1,
   %%%%%%%%%%%%%%%%
   proteina(ID,Primary,Secondary,Reference_aa_input_prot),
   write(Primary),nl,
   convert_list(Primary,Codes),
   write(Codes),nl,
   statistics(runtime,_),
   (options(mutation(List)),
     mutate(Primary,Primary1,List),
     retractall(new_protein(_,_,_,_)),
     assert(new_protein(new,Primary1,Secondary,Reference_aa_input_prot))
    ;
     \+options(mutation(_List)),
     Primary1=Primary
    )
    ,
    (options(use_original_tertiary(_)),  %% if needed original positions -> get them (in the correct current model)
      prepare_configuration(originalTertiary,Primary,Secondary,CommandsCodeT,CommandsTertiaryT),
      originalCodes(RequestedCodes,_),
      constrain(Primary,[],TertiaryT,CodeT,PrefDomT,CF,_PlacedTemplates,OriginalMatrices,
              _TotalEnergy,_TorsionEnergy,ContactEnergyT,DistanceListT,ContactEnergyStructuredT,
              0,_TorsContributions,_CorrContributions,
              _OriTertiary,_,RequestedCodes,CommandsCodeT,CommandsTertiaryT),
      my_labeling(CodeT,TertiaryT,_PCost,_Src,PrefDomT,_TotalEnergy,_TorsionEnergy,ContactEnergyT,DistanceListT,ContactEnergyStructuredT),
      TertiaryT=[OriginalTertiaryCA,_],
      %% now load the mutation db (if working with mutated sequence)
     (options(mutation(_)),!,
      init_data(new);
      true
     )
    ;
    \+options(use_original_tertiary(_))),

   prepare_configuration(noOriginalTertiary,Primary1,Secondary,CommandsCode,CommandsTertiary),
   originalCodes(OriginalCodes,_),
   constrain(Primary1,Secondary,Tertiary,Code,PrefDom,CF,PlacedTemplates,PlacedMatrices,
             TotalEnergy,TorsionEnergy,ContactEnergy,DistanceList,ContactEnergyStructured,
             InputEnergy,TorsContributions,CorrContributions,
             OriginalTertiaryCA,OriginalMatrices,OriginalCodes,CommandsCode,CommandsTertiary),
   !,
   write(constrainok),nl,
   constraint_time,
   %%%%%%%%%%%%%%%%%
   my_labeling(Code,Tertiary,PCost,Src,PrefDom,TotalEnergy,TorsionEnergy,ContactEnergy,DistanceList,ContactEnergyStructured),
   %%%%%%%% NEW: PATCH TO INCLUDE THE ORIENTATION ENERGY  February 10, 2011.
   Tertiary=[CAS,_],
   orientation_energy(CAS,ORE), OrientationEnergy is integer(ORE),
   TOT is OrientationEnergy+TotalEnergy,
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   format("### Torsion Energy: ~d,  Contact Energy: ~d, Orientation Energy: ~d~n### Total Energy: ~d\n",
          [TorsionEnergy,ContactEnergy,OrientationEnergy,TOT]),
   search_time(Time).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% add_sequence asserts the primary sequence of a protein
%%%    temporarily addressed as "new"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%10/02/2011

add_sequence(Primary) :-
   retractall(new_protein(_,_,_,_)),
   assert(new_protein(new,Primary,[],1)).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%       CONSTRAIN PHASE      %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

constrain(Primary,Secondary, Tertiary, Code, PrefDom,CF,PlacedTemplates,PlacedMatrices,
          TotalEnergy,TorsionEnergy,ContactEnergy,DistanceListStructured,
          ContactEnergyStructured,InputEnergy,TorsContributionsS,CorrContributionsS,
          OriginalTertiary,OriginalMatrices,OriginalCodes,CommandsCode,CommandsTertiary) :-

   %%% Retrieve parameters
   datalength(DL), identity_int(MatID),
   resolution(Res),
   min_distance(MD), min_distance_cg(MDcg),
   length(Primary,N), P is 3*N,
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%% Min distance^2 between non consecutive Calphas
   MDI_2 is integer(MD*Res*MD*Res),
   MDIcg_2 is integer(MDcg*Res*MDcg*Res),
   %%% diameter^2
   DIAM_2 is integer(CF*Res*CF*Res), %%% diameter^2
   %%% Split of the tertiary lists
   Tertiary=[CA,Centroids],

   %%% Generate FD variables
   length(CA,P),         %%% X,Y,Z where each tuple begin
   length(Centroids,P),  %%% X,Y,Z where each tuple begin

   %%% Assign Domains
   write('Codes restrictions:'),nl,
   domain_tuple(DL,Primary,Code,PrefDom, [First|RestTemplist],1,CommandsCode,OriginalCodes),

   ( %%%% OPTIONS - if then
     options(box(AACode, [RefX,RefY,RefZ], [[MinX,MaxX],[MinY,MaxY],[MinZ,MaxZ]])),!,
     box_constraint(CA, AACode, [RefX,RefY,RefZ], [[MinX,MaxX],[MinY,MaxY],[MinZ,MaxZ]])
     ; true ),

     ( %%%% OPTIONS - if then
     options(ellipsoid(AACode, [RefX,RefY,RefZ], [[F1x,F1y,F1z],[F2x,F2y,F2z],Rad])),!,
     ellipsoid_constraint(CA, AACode, [RefX,RefY,RefZ], [[F1x,F1y,F1z],[F2x,F2y,F2z],Rad])
     ; true ),



   First=[FirstT,[_C1x,_C1y,_C1z,C2x,C2y,C2z,C3x,C3y,C3z,_C4x,_C4y,_C4z]],
   prefix(CA,FirstT), %%% The first template is a prefix of Backbone (CA)
   Centroids=[_X1,_Y1,_Z1|RestCentroids],
   prefix(RestCentroids,[C2x,C2y,C2z,C3x,C3y,C3z]), %%% The first template is a prefix of Centroids (place only 2 central centroids)
   next_template(Code,MatrixList),  %%% Global constraint: only compatible consecutive templates
                                    %%% Matrix list is the list of rotation matrix that may link one template to the next one
   %%%% DISTANCE PART (on CA)
   distance(CA,DistanceListStructured,1,Secondary), %% also impose constraint on ssbonds
   split_matrix(DistanceListStructured,FlatALL,_NEXTs),
   %%% all_distant:
   domain(FlatALL,MDI_2,DIAM_2),

   ( %%%% OPTIONS - if then
    options(use_centroids),!,
    distance(Centroids,DistanceListStructuredCentroids,1,[]), %% distance between centroids
    split_matrix(DistanceListStructuredCentroids,FlatALLC,_),
    domain(FlatALLC,MDIcg_2,DIAM_2),
    distance_ca_centroids(CA,Centroids,DistanceListStructuredCACentroids,1), %% distance between CA and centroids
    split_matrix(DistanceListStructuredCACentroids,FlatALLCAC,_),
    domain(FlatALLCAC,MDI_2,DIAM_2)
    ;  true ),

   %%% CONTIGUOUS CONSTRAINTS. Actually this part is commented
   %%% Min and Max allowed distance^2 between consecutive Calphas
   %next_CaCadistance(NextD),
   %NEXT_min_2 is integer(NextD*Res*NextD*Res*0.90),
   %NEXT_MAX_2 is integer(NextD*Res*NextD*Res*1.10),
   %domain(NEXTs,NEXT_min_2,NEXT_MAX_2),
   %triangular_disequality(DistanceListStructured),
   %structural_distance(DistanceListStructured),

   %%% propagate position of templates and rotation matrices
   %%% based on commands CommandsTertiary (if requested to use original position, handle correctly position and orientation)
   link_templates(Tertiary,CommandsTertiary,MatID,RestTemplist,MatrixList,
           RestPlacedTemplates, PlacedMatrices,OriginalMatrices,OriginalTertiary),
   PlacedTemplates=[First|RestPlacedTemplates],

   %%%% ENERGY PART
   TotalEnergy #< InputEnergy,
   torsional_energy(Primary,Code, TorsContributions,CorrContributions),
   append(TorsContributions,CorrContributions,Encontribs),
   sum(TorsContributions,#=,TorsContributionsS),
   sum(CorrContributions,#=,CorrContributionsS),
   sum(Encontribs,#=,TorsionEnergy),

   (%%%% OPTIONS - if then else
     options(use_centroids),!,
     %% contribution from CG - CG contact
     contact_energy_centr(Primary,DistanceListStructuredCentroids,ContactEnergy1,ContactEnergyStructured),
     %% contribution from CA - CA contact
     contact_energy_ca(DistanceListStructured,ContactEnergy2,_ContactEnergyStructured1),
 %    ContactEnergy#=ContactEnergy1+ContactEnergy2
  %%%% NEW: weighted sum
     ContactEnergy#=((ContactEnergy1+ContactEnergy2)*4)/10
     ; %%ELSE: no centroids
     contact_energy(Primary,DistanceListStructured,ContactEnergy,ContactEnergyStructured)),

   TotalEnergy#=TorsionEnergy+ContactEnergy,
   !.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  domain_tuple sets the admissible "code"
%%%  value for all the variables representing tuples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 14/11/2009
%%% if the template is in the DB it computes the preferred
%%% ordered domains (from most probable to least probable)
%%% otherwise it assign a pseudo template [-1,-1,-1,...]
%%% and do the same.
%%% Moreover, admissible values for corresponding
%%% tuple variables are assigned in a table combinatorial constraint
%%% NEW: if the code is included in a secondary structure element, force it to the corresponding template ID
%%% NEW: Output T is a pair of constrained variables [Templ,Centroids]
%%%      where Templ is a list of 3d coordinates of CAlpha and Centroids the corresponding list of centroids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 24/11/2009

domain_tuple(DATALENGTH,PRIMARY,[],[],[],_,_CommandsCode,_):-
      length(PRIMARY,N),
      N < DATALENGTH,!,nl.
domain_tuple(DATALENGTH,[P|RIMARY],[C|Code],[ListDom|RDom],[T|RT],N,[Command|CommandsCode],OriginalCodes) :-
write(Command),
      length(PREF,DATALENGTH),
      append(PREF,_,[P|RIMARY]),
      %% test if N in a helix
      (
       Command='h',!,
       cleantupla([-2,-2,-2,-2],CA,[Centroids2,Centroids3],_,Alpha,_),
       PREF=[_,AA2,AA3,_],
       member([AA2,CPoint2],Centroids2),  %% select actual amino acid
       member([AA3,CPoint3],Centroids3),
       CA=[X1,Y1,Z1,_,_,_,_,_,_,X4,Y4,Z4],
       append([X1,Y1,Z1],CPoint2,CentroidsT1),
       append(CentroidsT1,CPoint3,CentroidsT2),
       append(CentroidsT2,[X4,Y4,Z4],Centroids),
       T=[CA,Centroids],
       C=Alpha,
       ListDom=[Alpha]
       ;
      %% test if N in a strand
       Command='s',!,
       cleantupla([-3,-3,-3,-3],CA,[Centroids2,Centroids3],_,Beta,_),
       PREF=[_,AA2,AA3,_],
       member([AA2,CPoint2],Centroids2),  %% select actual amino acid
       member([AA3,CPoint3],Centroids3),
       CA=[X1,Y1,Z1,_,_,_,_,_,_,X4,Y4,Z4],
       append([X1,Y1,Z1],CPoint2,CentroidsT1),
       append(CentroidsT1,CPoint3,CentroidsT2),
       append(CentroidsT2,[X4,Y4,Z4],Centroids),
       T=[CA,Centroids],
       C=Beta,
       ListDom=[Beta]
      ;
      %% normal
      TreD is DATALENGTH*3,
      convert_list(PREF,CL),
      (
         findall(Found,
         (cleantupla(CL,_,_,_,Id,_),
          (ground(OriginalCodes),

             (Command\='c',                  %% if not original code needed, check that original codes are not in the way...
              \+member(Id,OriginalCodes);
             true)
          ;
          true),
          Found=ok
         ),
         [_|_]),
         !, %%% known tuple (and not in the original code list
         CODE=CL;
         length(CODE,DATALENGTH),   %%% unknown tuple
         same_value(CODE,-1)        %%% special case
      ),
      findall(Id, cleantupla(CODE,_,_,_Freq,Id,_), ListTemp),

      %% remove from the list the original codes for the protein (avoid using templates from the same protein)
      (Command='c',
       nth1(N,OriginalCodes,OriC),
       list_to_fdset([OriC],D);  %% only original code
       Command='-',
       list_to_fdset(ListTemp,DTemp),

       (ground(OriginalCodes),!,              %%% no info about original codes -> dummy domain
        list_to_fdset(OriginalCodes,NoTemp);  %% remove original code from the list (cannot use original fragment)
        list_to_fdset([-1],NoTemp)
       ),
       fdset_subtract(DTemp,NoTemp,D)
      ),
%write(a(D,CL,CODE,ListTemp)),nl,
      C in_set D,
      findall(Freq-Id, cleantupla(CODE,_,_,Freq,Id,_), List1),
      rev_sort(List1, RSList1),
      findall(B, member(_A-B,RSList1), ListDom),
      %%% part from create_templ:
      findall([ID|Templ], cleantupla(CODE,Templ,_Centroids,_,ID,_), PossibleTemplates),
      length(T1,TreD),
      table([[C|T1]],PossibleTemplates),
      findall([ID|Centroids1], (
          cleantupla(CODE,CA,[Centroids2,Centroids3],_,ID,_),
          PREF=[_,AA2,AA3,_],
          member([AA2,CPoint2],Centroids2),  %% select actual amino acid
          member([AA3,CPoint3],Centroids3),
          CA=[X1,Y1,Z1,_,_,_,_,_,_,X4,Y4,Z4],
          append([X1,Y1,Z1],CPoint2,CentroidsT1),
          append(CentroidsT1,CPoint3,CentroidsT2),
          append(CentroidsT2,[X4,Y4,Z4],Centroids1)
      ), PossibleTemplatesC),
      length(T2,TreD),
      table([[C|T2]],PossibleTemplatesC),
      T=[T1,T2]
      ),
%write(a),nl,
      N1 is N+1,
      domain_tuple(DATALENGTH,RIMARY,Code,RDom,RT,N1,CommandsCode,OriginalCodes).



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% next_template uses "next" to force arc consistency
%%%% using combinatorial constraint "table"
%%%% between consecutive tuplets
%%%% Moreover, a table of possible matrix rotation
%%%% is stored.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%14/11/2009

next_template([],[]).
next_template([_],[]).
next_template([C1,C2|R],[T|MList]) :-
   fd_set(C1,Set1),  fdset_to_list(Set1,List1),
   fd_set(C2,Set2),  fdset_to_list(Set2,List2),
   create_pairs(List1,List2,TABLE),
   table([[C1,C2]],TABLE), %%% SICStus BUILT-IN
   %%%% cut from create_matrix
   findall([LC1,LC2|M],
       (member(LC1,List1), member(LC2,List2), next(LC1,LC2,M)),
            PossibleMatrices),
   length(T,9), %%% Rotation matrix (3x3)
   table([[C1,C2|T]],PossibleMatrices),
   next_template([C2|R],MList).

create_pairs([],_,[]).
create_pairs([A|R],List,TABLE) :-
   create_pairs(R,List,T),
   findall([A,B],(cleantupla(_,_,_,_,B,_),member(B,List),next(A,B,_)),TABLE,T).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare_configuration/5 :
%%% given the flag HandleTertiary (force to use original codes as preprocessed in the db
%%% originalTertiary = use original codes from db to generate the corresponding tertiary model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prepare_configuration(HandleTertiary,Primary,Secondary,Commands,CommandsTertiaryCopy):-
   length(Primary,N),
   N3 is N-3,
   length(Commands0,N3),
   length(CommandsTertiaryCopy0,N),

   %%% requested to fix all codes -> original protein production
   (
    (HandleTertiary=originalTertiary;
     options(use_original_codes)),
    reference_aa_input_prot(Ofs),
    Start is 1+Ofs-1,
    End is N+Ofs-1,
    updateCommands([[Start,End]],Commands0,Commands1,c,codes)
    ;
    HandleTertiary\=originalTertiary,
    \+options(use_original_codes),
    Commands1=Commands0
   ),

   %%
   (options(use_original_codes(Ranges1)),
    updateCommands(Ranges1,Commands1,Commands21,c,codes)
    ;
    \+options(use_original_codes(_Ranges)),
    Commands1=Commands21
   ),

   (HandleTertiary\=originalTertiary,
    options(use_original_tertiary(Ranges0)),  %% not in the pre-calculation of tertiary -> based on tertiary -> also need original codes!
    updateCommands(Ranges0,Commands21,Commands2,c,codes)
    ;
    (HandleTertiary=originalTertiary;  %% no need to equate CAs, computes all now
     \+options(use_original_tertiary(_Ranges))),
     Commands21=Commands2
    ),


   (options(use_secondary),
    findall([S,E],member(strand(S,E),Secondary),StrandRanges),
    updateCommands(StrandRanges,Commands2,Commands3,s,codes),
    findall([S,E],member(helix(S,E),Secondary),HelixRanges),
    updateCommands(HelixRanges,Commands3,Commands4,h,codes)
    ;
    \+options(use_secondary),
     Commands2=Commands4
    ),

   reference_aa_input_prot(Ofs),
    Start is 1+Ofs-1,
    End is N+Ofs-1,
    updateCommands([[Start,End]],Commands4,Commands5,'-',codes),

   (HandleTertiary\=originalTertiary,
    options(use_original_tertiary(Ranges0)),  %% not in the pre-calculation of tertiary -> copy the needed tertiary
    updateCommands(Ranges0,CommandsTertiaryCopy0,CommandsTertiaryCopy1,t,tertiary)
    %%% simple copy from ranges: here referred to amino acids (not codes)
    %%%withdraw_original_tertiary(...)
    ;
    (HandleTertiary=originalTertiary;  %% no need to equate CAs, computes all now
     \+options(use_original_tertiary(_Ranges))),
     CommandsTertiaryCopy0=CommandsTertiaryCopy1
     ),
    updateCommands([[Start,End]],CommandsTertiaryCopy1,CommandsTertiaryCopy,'-',tertiary),

   %% tertiary copy -> referred to codes (when to link tertiary or detatch with an approximate box landing

   Commands5=Commands,
   write(Commands),nl,
   write(CommandsTertiaryCopy),nl.

%% patches the variables: if already a value is present, keeps it!
%% this way preserve the priority of commands
%% start, end: refers to amino acids in the original reference

updateCommands([],Commands,Commands,_FillCommand,_).
updateCommands([[Start,End]|Rest],Commands,Commands1,FillCommand,CodesOrTertiary):-
  reference_aa_input_prot(Ofs),
  Pos1 is Start-Ofs+1,
  (CodesOrTertiary=codes,
   Pos2 is End-Ofs+1 -3;  %% from aa pos to code pos (must be reduced by 3)
   CodesOrTertiary=tertiary,
   Pos2 is End-Ofs+1),    %% from aa pos to tertiary pos
  Pos1m1 is Pos1-1,

  Len is Pos2-Pos1+1,
  (Len>0,!,
   length(Patch,Len),
   same_value(Patch,FillCommand),

  sublist(Commands,Prefix,0,Pos1m1,_),
  sublist(Commands,Previous,Pos1m1,Len,_),
  sublist(Commands,Suff,Pos2,_,0),

  patch(Previous,Patch,Patched),

  append(Prefix,Patched,T1),
  append(T1,Suff,Commands2)
  ;
  Len<1, Commands2=Commands),
  updateCommands(Rest,Commands2,Commands1,FillCommand,CodesOrTertiary).



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% distance/4 retrieves from Tertiary = [X1,Y1,Z1,...,XN,YN,ZN]
%%%% a diagonal distance matrix done of lists like this:
%%%% [i,[D_(i,i+1),...,D_(i,n)]] for i = 1, ..., N-2
%%%% where D_(i,j) is the squared Euclidean distance betweeen
%%%% Xi,Yi,Zi and Xj,Yj,Zj.
%%%% In the same time, it adds distance constraints between
%%%% the known ssbonds stored in "Sec".
%%%% Aux: distance/8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%14/11/2009

distance([],[],_,_):-!.
distance([_,_,_],[],_,_):-!.
distance([Xi,Yi,Zi|Tertiary],[[I, ROW ]|TriangDist],I,Sec):-
     I1 is I+1,
     distance(Tertiary, TriangDist, I1,Sec),
     distance(Xi,Yi,Zi,Tertiary,ROW ,I, I1,Sec).

distance(_,_,_,[],[],_,_,_).
distance(Xi,Yi,Zi,[Xj,Yj,Zj|Tertiary],[D|Row],I,J,Sec):-
    D #= (Xi-Xj)*(Xi-Xj)+(Yi-Yj)*(Yi-Yj)+(Zi-Zj)*(Zi-Zj),
    (
      reference_aa_input_prot(Ofs),
      I1 is I-1+Ofs,
      J1 is J-1+Ofs,
      member(ssbond(I1,J1),Sec),!,
     resolution(Res),
     distanceSsbond(SSD),
     SsbondDist  is integer(Res*SSD),
     SsbondDist_2 is integer(Res*SSD*Res*SSD),
     D #=< SsbondDist_2,
     abs(Xi-Xj) #=< SsbondDist,
     abs(Yi-Yj) #=< SsbondDist,
     abs(Zi-Zj) #=< SsbondDist,
     format("### added ssbond constraint ~d-~d\n",[I,J]);
      true),
     J2 is J+1,
     distance(Xi,Yi,Zi,Tertiary,Row,I,J2,Sec).

distance_ca_centroids([],_,[],_):-!.
distance_ca_centroids([_,_,_],_,[],_):-!.
distance_ca_centroids([Xi,Yi,Zi|Tertiary],Centroids,[[I, ROW ]|TriangDist],I):-
     I1 is I+1,
     distance_ca_centroids(Tertiary, Centroids,TriangDist, I1),
     distance_ca_centroids(Xi,Yi,Zi,Centroids,ROW ,I, 1).

distance_ca_centroids(_,_,_,[],[],_,_).
distance_ca_centroids(Xi,Yi,Zi,[Xj,Yj,Zj|Centroids],[D|Row],I,J):-
    (I#\=J,!,  %% skip check on the same aa ca and centroid!)
    D #= (Xi-Xj)*(Xi-Xj)+(Yi-Yj)*(Yi-Yj)+(Zi-Zj)*(Zi-Zj);
    true),
    J1 is J+1,
    distance_ca_centroids(Xi,Yi,Zi,Centroids,Row,I,J1).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% triangular_disequality/1 (aux triangular_disequality_row/3)
%%%% It runs on the triangular matrix DLS to call
%%%% the auxiliary predicates row_search.
%%%% row_search adds a triangular constraint.
%%%% Given three points I, J, K, then
%%%%   D_IJ  #=< 2*(D_IK + D_JK),
%%%%   D_IK  #=< 2*(D_IJ + D_JK),
%%%%   D_JK  #=< 2*(D_IK + D_IJ);
%%%% Observe that they are squared.
%%%% We know that sqrt(a) =< sqrt(b) + sqrt(c)
%%%% Therefore a =< b + c + 2 sqrt(b) sqrt(c)
%%%% We set a =< 2 (b + c).
%%%% This constraint is an approximation that
%%%% however allows propagation.
%%%% The problem is that it adds a number of constraints
%%%% cubic wrt protein length !
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 14/11/2009

triangular_disequality([]).
triangular_disequality([[_,_]]):-!.
triangular_disequality([[I,Row]| DLS]):-
    I1 is I+1,
    triangular_disequality_row(I1,Row,DLS),
    triangular_disequality(DLS).

triangular_disequality_row(_J,[],_DLS).
triangular_disequality_row(_J,[_],_DLS).
triangular_disequality_row( J,[D_IJ,D_IJ1|R],DLS):-
    J1 is J+1,
    triangular_disequality_row(J,D_IJ,J1,[D_IJ1|R],DLS),
    triangular_disequality_row(J1,[D_IJ1|R],DLS).

triangular_disequality_row(_,_,_,[],_).
triangular_disequality_row(J,D_IJ,J1,[D_IJ1|Row],DLS):-
   row_search(J,D_IJ,J1,D_IJ1,DLS),
   J2 is J1+1,
   triangular_disequality_row(J,D_IJ,J2,Row,DLS).

%%%% Note: point "I" is implicit

row_search(J,D_IJ, K,D_IK, DLS):-
   (member([J,ListJ],DLS),
    JK is K-J, nth1(JK,ListJ,D_JK),
    !,
    D_IJ  #=< 2*(D_IK + D_JK),
    D_IK  #=< 2*(D_IJ + D_JK),
    D_JK  #=< 2*(D_IK + D_IJ);
   true). %%% AGO: dubbio. Ma quand'e' che puo' fallire?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% structural_distance/1 (aux structural_distance/3)
%%%% It runs on the triangular matrix DLS
%%%% Every pair I,J and squared distance D_IJ  is addressed.
%%%% Added the constraint
%%%% D_IJ #< min(Diameter, 3.8*abs(I-J))^2
%%%% using the correct scale/resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 14/11/2009

structural_distance([]).
structural_distance([[I,ROW]|DLS]):-
  I1 is I+1,
  structural_distance_row(I,I1,ROW),
  structural_distance(DLS).

structural_distance_row(_ ,_,[]).
structural_distance_row(I , J ,[D_IJ|ROW]):-
   resolution(Res),
   distanceCaCa(CaDist),
   SD is integer(CaDist*Res),
   distanceCaCaStep2(CaDist2),
   SD2 is integer(CaDist2*Res),
   Delta is J-I,
   Delta2 is Delta//2,
   DeltaM is Delta rem 2,
   Dist is Delta2*SD2+DeltaM*SD,
   max_distance(Max),
   (Dist < Max*Res,!,
    D_IJ #=< Dist*Dist
    ;
    true),
   J1 is J+1,
   structural_distance_row(I, J1,ROW).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% box_constraint/4
%%%%  retrieves the Offset
%%%%  and forces the coordinated of the AACode-Offset
%%%%  to stay inside the box
%%%%  identified by the vertexes of the 4th parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%08/02/2010

box_constraint(CA,AACode, [RefX,RefY,RefZ],
         [[MinX,MaxX],[MinY,MaxY],[MinZ,MaxZ]]):-
    %%% write(setBox),nl,
    reference_aa_input_prot(Ofs),
    Posx is 3*(AACode-Ofs), Posy is 3*(AACode-Ofs)+1,
    Posz is 3*(AACode-Ofs)+2,
    resolution(Res),
    nth0(Posx,CA,AAx), nth0(Posy,CA,AAy), nth0(Posz,CA,AAz),
    Xm is integer(Res*(MinX-RefX)), Ym is integer(Res*(MinY-RefY)),
    Zm is integer(Res*(MinZ-RefZ)),
    XM is integer(Res*(MaxX-RefX)), YM is integer(Res*(MaxY-RefY)),
    ZM is integer(Res*(MaxZ-RefZ)),
    %%% Box
    AAx #>= Xm, AAy #>= Ym, AAz #>= Zm,
    AAx #=< XM, AAy #=< YM, AAz #=< ZM.

ellipsoid_constraint(CA,AACode, [RefX,RefY,RefZ],
         [[F1x,F1y,F1z],[F2x,F2y,F2z],Rad]):-
    %%% write(setBox),nl,
    reference_aa_input_prot(Ofs),
    Posx is 3*(AACode-Ofs), Posy is 3*(AACode-Ofs)+1,
    Posz is 3*(AACode-Ofs)+2,
    resolution(Res),
    nth0(Posx,CA,AAx), nth0(Posy,CA,AAy), nth0(Posz,CA,AAz),


    F1xI is integer(Res*(F1x-RefX)),
    F1yI is integer(Res*(F1y-RefY)),
    F1zI is integer(Res*(F1z-RefZ)),
    F2xI is integer(Res*(F2x-RefX)),
    F2yI is integer(Res*(F2y-RefY)),
    F2zI is integer(Res*(F2z-RefZ)),

    A#=(AAx-F1xI)*(AAx-F1xI)+
    (AAy-F1yI)*(AAy-F1yI)+
    (AAz-F1zI)*(AAz-F1zI)+
    (AAx-F2xI)*(AAx-F2xI)+
    (AAy-F2yI)*(AAy-F2yI)+
    (AAz-F2zI)*(AAz-F2zI),


    %%% Ellipsoid
    A #=< Rad*Rad*Res*Res.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% link_templates/5 generates the constraints
%% that relates the 3 last points of a tuple
%% with the 3 first points of the consecutive
%% one as well as the fourth point of the second
%% as effect of distance minimiaztion and rotation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 16/11/2009
%%  PlTempl is a list of rotated and shifted templates,
%%      used for debugging
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

link_templates(_,_,_OldMatrix,[],[],[],[],[],_).
link_templates([
         [X1,Y1,Z1,    X2,Y2,Z2,    X3,Y3,Z3,    X4,Y4,Z4,    X5,Y5,Z5    |CATertiary],
         [_X1C,_Y1C,_Z1C, X2C,Y2C,Z2C, X3C,Y3C,Z3C, X4C,Y4C,Z4C, X5C,Y5C,Z5C |CentrTertiary]
         ],
         [T1,T2,T3,T4,T5|CommandsTer],
    OldMatrix,
    [[Templ,TemplC]|RestTempl],[Matrix|MatrixList],[[PlTempl,PlTemplC]|RestPlTempl],
    [NewMatrix|RestMatrices],[OriginalMatrix|RestOriginalMatrices],
    [Ox1,Oy1,Oz1,Ox2,Oy2,Oz2,Ox3,Oy3,Oz3,Ox4,Oy4,Oz4,Ox5,Oy5,Oz5|RestOriginalPos]):-

    (T1='-',T2='t',!,             %% in case of joining a copied pattern, link original reference
     NewMatrix=OriginalMatrix
     ;
     matrix_mult_int(Matrix,OldMatrix,NewMatrix) %% otherwise, simply connect with current next rotation matrix
    ),

    rotate_coord_int(Templ, NewMatrix, [VX2,VY2,VZ2,    VX3,VY3,VZ3,    VX4,VY4,VZ4,    VX5,VY5,VZ5]),
    rotate_coord_int(TemplC,NewMatrix, [VX2C,VY2C,VZ2C, VX3C,VY3C,VZ3C, VX4C,VY4C,VZ4C, VX5C,VY5C,VZ5C]),
    Dx #= X4-VX4, Dy #= Y4-VY4, Dz #= Z4-VZ4,

    shift_int([VX5,VY5,VZ5], Dx, Dy, Dz ,[VX5Candidate,VY5Candidate,VZ5Candidate]),
    resolution(Res),
    toleranceLinkDist(TD),
    ToleranceDist2 is integer(Res*TD*2),
    ToleranceDist15 is integer(Res*TD*1.5),
    % ToleranceDist1 is integer(Res*TD),
    (T4='-',T5='t',!,
     ToleranceDist15 #> abs(VX5Candidate-X5),
     ToleranceDist15 #> abs(VY5Candidate-Y5),
     ToleranceDist15 #> abs(VZ5Candidate-Z5)
     ;
     T3='-',T4='t',!,
%    true
     ToleranceDist2 #> abs(VX5Candidate-X5),
     ToleranceDist2 #> abs(VY5Candidate-Y5),
     ToleranceDist2 #> abs(VZ5Candidate-Z5)
     ;
     T2='-',T3='t',!,
     true
%    ToleranceDist2 #> abs(VX5Candidate-X5),
%    ToleranceDist2 #> abs(VY5Candidate-Y5),
%    ToleranceDist2 #> abs(VZ5Candidate-Z5)
     ;
     T1='-',T2='t',!,
     true
%    ToleranceDist2 #> abs(VX5Candidate-X5),
%    ToleranceDist2 #> abs(VY5Candidate-Y5),
%    ToleranceDist2 #> abs(VZ5Candidate-Z5)
     ;
     X5=VX5Candidate,
     Y5=VY5Candidate,
     Z5=VZ5Candidate),


    shift_int([VX4C,VY4C,VZ4C], Dx, Dy, Dz ,[X4C,Y4C,Z4C]),
    shift_int([VX2,VY2,VZ2, VX3,VY3,VZ3, VX4,VY4,VZ4, VX5,VY5,VZ5],
               Dx, Dy, Dz ,PlTempl),
    shift_int([VX2C,VY2C,VZ2C, VX3C,VY3C,VZ3C, VX4C,VY4C,VZ4C, VX5C,VY5C,VZ5C],
               Dx, Dy, Dz ,PlTemplC),

     %% map original positions and ranges
    (T1='t',!,
     X1=Ox1, Y1=Oy1, Z1=Oz1
     ;
     T2='t',!,
     X2=Ox2, Y2=Oy2, Z2=Oz2
     ;
     true),


    link_templates([
    [X2,Y2,Z2,    X3,Y3,Z3,    X4,Y4,Z4,    X5,Y5,Z5    |CATertiary],
    [X2C,Y2C,Z2C, X3C,Y3C,Z3C, X4C,Y4C,Z4C, X5C,Y5C,Z5C |CentrTertiary]
    ],
    [T2,T3,T4,T5|CommandsTer],
    NewMatrix,RestTempl,MatrixList,RestPlTempl,RestMatrices,RestOriginalMatrices,
    [Ox2,Oy2,Oz2,Ox3,Oy3,Oz3,Ox4,Oy4,Oz4,Ox5,Oy5,Oz5|RestOriginalPos]).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   ENERGY CONSTRAINTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% contact_energy/4 (aux: contact_energy_row/4)
%%% Given the primary sequence and the triangular
%%%   distance matrix returns the constrained contact energy, and
%%%   the energy row by row.
%%% If the option is active the contact between Centroids is
%%% also computed
%%% It calls tablep, the 20x20 potential matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 08/02/2010

%% model with only CA
contact_energy(_,[],0,[]).
contact_energy([AA|Rest],[[N,L]|DistanceListStructured],ContactEnergy,
   [[N,En]|ContactEnergyStr]):-
   ContactEnergy #= ContactEnergy1+ContactEnergy2,
   contact_energy_row(AA,Rest,L,ContactEnergy2,En), %%% AGO DUBBIO su Min in questo caso
   contact_energy(Rest,DistanceListStructured,ContactEnergy1,ContactEnergyStr).

%% model with centroids and ca
%% skip first and last centroid (by construction are meaningless)
contact_energy_centr(_,[],0,[]).
contact_energy_centr([AA|Rest],[[N,L]|DistanceListStructured],ContactEnergy,[[N,En]|ContactEnergyStr]):-
   (N>1 ->
     (ContactEnergy #= ContactEnergy1+ContactEnergy2,
      %%remove last centroid (and corr. aa) because is meaningless
      append(Rest1,[_],Rest), append(L1,[_],L),
      contact_energy_row(AA,Rest1,L1,ContactEnergy2,En));
    N=1 -> ContactEnergy #= ContactEnergy1),
   contact_energy_centr(Rest,DistanceListStructured,ContactEnergy1,ContactEnergyStr).

contact_energy_ca([],0,[]).
contact_energy_ca([[N,L]|DistanceListStructured],ContactEnergy,[[N,En]|ContactEnergyStr]):-
   ContactEnergy #= ContactEnergy1+ContactEnergy2,
   contact_energy_row_ca(L,ContactEnergy2,En),
   contact_energy_ca(DistanceListStructured,ContactEnergy1,ContactEnergyStr).

%% contribution as if it were an asn aminoacid

contact_energy_row_ca([],0,[]).
contact_energy_row_ca([Dist|RestDist],ContactEnergy1,[Contrib|RestContr]):-
   tablep(n, n, Pot),
   resolution(Res),
   Min is 4.8,         %% cutoff distance for CA-CA pair
   MinDist is integer(Min*Res*Min*Res),
   ConstantContrib is integer(Pot * MinDist / 100),
   (Dist #> MinDist) #=>
      (Dist100 #= Dist/100 #/\  Contrib #= ConstantContrib/Dist100),
  (Dist #=< MinDist ) #=> (Contrib #=  Pot),
  ContactEnergy1 #= Contrib + ContactEnergy2,
  contact_energy_row_ca(RestDist,ContactEnergy2,RestContr).


contact_energy_row(_,_,[],0,[]).
contact_energy_row(AA,[AAJ|Rest],[Dist|RestDist],ContactEnergy1,[Contrib|RestContr]):-
   tablep(AA, AAJ, Pot),
   resolution(Res),
   %%% The radii of the centroids are computed and summed
   radius(AA,R1), radius(AAJ,R2),
   Min is R1+R2,
   MinDist is integer(Min*Res*Min*Res),
   ConstantContrib is integer(Pot * MinDist / 100),
   (Dist #> MinDist) #=>
     (Dist100 #= Dist/100 #/\ Contrib #= ConstantContrib/Dist100),
   (Dist #=< MinDist ) #=> (Contrib #=  Pot),
   ContactEnergy1 #= Contrib + ContactEnergy2,
   contact_energy_row(AA,Rest,RestDist,ContactEnergy2,RestContr).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% torsional_energy/4
%%%  it computes two sets of potential values
%%%  the first for the torsions
%%%  the second for the correlations between successive tuples
%%% it calls
%%%  - "prof" that returns the list of profiles for the two codes
%%%       C1 and C2 associated to B,C
%%%  - "corr" that returnd the correlation list
%%%  - "tors_templ".
%%% These facts are stored in the file  PROTNAMEcurrenttorsdb.pl
%%% generated during the preprocessing stage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%08/02/2010

torsional_energy(_,[], [],[]).
torsional_energy([_A,B,C,_D],[Code], [T],[]):-

%%write([B,C]),nl,
     %convert_list([B,C],[C1,C2]),
     B=C1,
     C=C2,
     fd_set(Code,Set1),  fdset_to_list(Set1,List1),
     prof(C1,C2,ListProf),
     findall([ID,E],(member(ID,List1), tors_templ(ID,T),
              member([T,E],ListProf)),
           TableList),
     table([[Code,T]],TableList).

torsional_energy([_A,B,C,D|R],[Code,Code1|R1], [T|TorsContr],[CorrContrib|Correlation]):-
     %convert_list([B,C],[C1,C2]),
%%write([B,C]),nl,
     C1=B,
     C2=C,
     fd_set(Code,Set1),  fdset_to_list(Set1,List1),
     prof(C1,C2,ListProf),
     findall([ID,E],(member(ID,List1), tors_templ(ID,T),
              member([T,E],ListProf)),
           TableList),
     table([[Code,T]],TableList),
     fd_set(Code1,Set2),  fdset_to_list(Set2,List2),
     corr(ListCorr),
     findall([ID1,ID2,E],
         (member(ID1,List1), member(ID2,List2),
          tors_templ(ID1,T1), tors_templ(ID2,T2),
          member([T1,T2,E],ListCorr)),
             CorrList),
     table([[Code,Code1,CorrContrib]],CorrList),
     torsional_energy([B,C,D|R],[Code1|R1],TorsContr,Correlation).


%%% FEBRUARY 2011
%%%%%%%%%%%% AGO. LA METTEREI COME "delta" alla fine.
%%%%%%%%%%%% SONO TROPPI I VINCOLI CHE INTRODUCE E CON TUTTE
%%%%%%%%%%%% QUELLE ROTAZIONI E RADICI NON CE NE FAREMMO MOLTO
%%%%%%%%%%%% PER QUESTO HO USATO "is" invece di #=
%%%%%%%%%%%% IN OGNI MODO PER USARE #= DOBBIAMO FAR IN MODO CHE TUTTO SIA INTERO...

orientation_energy(CAS,ORE) :-
   length(CAS,N3), M is N3//3-2, M1 is M+1,
   orientation_energy_aux(2,M,3,M1,CAS,ORE).

orientation_energy_aux(M,M,_,_,_,0) :-
   !.
orientation_energy_aux(I,M,N,N,CAS,ORE) :-
   !,
   I1 is I+1,
   J is I+2,
   orientation_energy_aux(I1,M,J,N,CAS,ORE).
orientation_energy_aux(I,M,J,N,CAS,ORE) :-
   J1 is J+1,
   orientation_energy_step(I,J,CAS,OR1),
   orientation_energy_aux(I,M,J1,N,CAS,OR2),
   ORE is OR1+OR2.

orientation_energy_step(I,J,CAS,Energy) :-
%%% FOGOLARI INSTRUCTIONS
%%% * Considero l'interazione fra CA[i] e Ca[j] dove i e j non sono ne' il
%%% primo ne' l'ultimo.
%%% * Calcolo la distanza fra CA[i] e CA[j] e il vettore distanza che
%%% chiamo v: v = CA[j] - CA[i]
%%% * Calcolo il vettore che va da CA[i-1] a CA[i] e quello che va
%%% da CA[i] a CA[i+1]. Li chiamo v1 e v2*/
%%%v1 = CA[i] - CA[i-1]
%%%v2 = CA[i+1] - CA[i]
%%% * Faccio lo stesso per CA[j]
%%% v5 = CA[j] - CA[j-1]
%%% v6 = CA[j+1] - CA[j]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute indices
   IZind is 3*I,    IYind is 3*I-1,  IXind is 3*I-2, %%% indices for CAI
   IZmind is 3*I-3, IYmind is 3*I-4, IXmind is 3*I-5, %%% indices for CAI-1
   IZpind is 3*I+3, IYpind is 3*I+2, IXpind is 3*I+1, %%% indices for CAI+1
   JZind is 3*J,    JYind is 3*J-1,  JXind is 3*J-2, %%% indices for CAJ
   JZmind is 3*J-3, JYmind is 3*J-4, JXmind is 3*J-5, %%% indices for CAJ-1
   JZpind is 3*J+3, JYpind is 3*J+2, JXpind is 3*J+1, %%% indices for CAJ+1
%%%% Retrieve coordinates
   nth1(IXind,CAS,XI),nth1(IYind,CAS,YI),nth1(IZind,CAS,ZI), %%% Coordinates for I
   nth1(IXmind,CAS,XIm),nth1(IYmind,CAS,YIm),nth1(IZmind,CAS,ZIm), %%% Coordinates for I-1
   nth1(IXpind,CAS,XIp),nth1(IYpind,CAS,YIp),nth1(IZpind,CAS,ZIp), %%% Coordinates for I+1
   nth1(JXind,CAS,XJ),nth1(JYind,CAS,YJ),nth1(JZind,CAS,ZJ), %%% Coordinates for J
   nth1(JXmind,CAS,XJm),nth1(JYmind,CAS,YJm),nth1(JZmind,CAS,ZJm), %%% Coordinates for J-1
   nth1(JXpind,CAS,XJp),nth1(JYpind,CAS,YJp),nth1(JZpind,CAS,ZJp), %%% Coordinates for J+1
%%%% Compute difference vectors
   VX  is XJ-XI, VY  is YJ-YI, VZ  is ZJ-ZI,       %% V vector = J - I
   V1X  is XI-XIm, V1Y  is YI-YIm, V1Z  is ZI-ZIm, %% V1 vector = I - (I-1)
   V2X  is XIp-XI, V2Y  is YIp-YI, V2Z  is ZIp-ZI, %% V2 vector = (I+1) - I
   V5X  is XJ-XJm, V5Y  is YJ-YJm, V5Z  is ZJ-ZJm, %% V5 vector = J - (J-1)
   V6X  is XJp-XJ, V6Y  is YJp-YJ, V6Z  is ZJp-ZJ, %% V2 vector = (J+1) - J
%%%%%%%%%%%%%%%%%%
%% * Calcolo il prodotto vettoriale fra v1 e v2 e lo chiamo v3 = v1 x v2
%% Questo e' un vettore ortogonale al piano su cui stanno v1 e v2 */
%% * Faccio lo stesso per v5 e v6 che ritorna v4 = v5 x v6
%%%%%%%%%%%%%%%
   prod_vect(V1X,V1Y,V1Z,V2X,V2Y,V2Z,V3X,V3Y,V3Z),
   prod_vect(V5X,V5Y,V5Z,V6X,V6Y,V6Z,V4X,V4Y,V4Z),
%%%%%%%%%%%%%%%%%%
%%  /* Normalizzo i vettori v, v3 e v4 */
%vnorm =   v/|v|
%v3norm = v3/|v3|
%v4norm = v4/|v4|
%%%%%%%%%%%%%%%%%%
  normalize(VX,VY,VZ,    VXnorm, VYnorm, VZnorm),
  normalize(V3X,V3Y,V3Z, VX3norm,VY3norm,VZ3norm),
  normalize(V4X,V4Y,V4Z, VX4norm,VY4norm,VZ4norm),
%%%%%%%%%%%%%%%%%%
%%% Sia d = distanza(CA[i], CA[j])
%%%%% NOTA AGO: se in clpfd useremmo S=d^2 per restare negli interi
%%%%%%%%%%%%%%%%%%
 % S #= VX*VX+VY*VY+VZ*VZ,
 % D*D #=< S, (D+1)*(D+1) #> S, %%% D e' la radice di S, piu' o meno
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 D is sqrt(VX*VX+VY*VY+VZ*VZ),

%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%
%%%/* Ora posso calcolare l'energia E */
%%% if d > 5.8
%%%    E  = 0
%%% else if (d > 5.6 && d < 5.8)
%%%    E  = -(5.8 - d)/(5.8 - 5.6)
%%% else if (d > 3.8 && d < 5.6)
%%%    E  = -1.0
%%%%%%%%% NOTA AGO. Assumo resolution 100.

  ( D > 580 ->  E = 0 ;
   (D > 560, D =< 580) -> E is (D-580)/20 ;
   D =< 560 -> E = -1 ),

%%%%%%%%%%%%%%
%%% /*Se E e' diverso da 0 procedo a scalarlo in base all'orientazione */
%%%
%%% if E != 0
%%% /* Calcolo il valore assoluto dei prodotti scalari v*v3, v*v4, v3*v4 */
%%% t = abs( (v*v3)*(v*v4)*(v3*v4) )
%%% if(t>0.7)
%%% E = E
%%% else if(t>0.5 && t < 0.7)
%%% E = E * (t - 0.5)/(0.7 - 0.5)
%%% else if (t<0.5)
%%% E = 0.0
%%% Il tutto prob. moltiplicato per 1000 visto i nostri parametri.

  prod_scal(VXnorm,VYnorm,VZnorm,VX3norm,VY3norm,VZ3norm,A1),
  prod_scal(VXnorm,VYnorm,VZnorm,VX4norm,VY4norm,VZ4norm,A2),
  prod_scal(VX3norm,VY3norm,VZ3norm,VX4norm,VY4norm,VZ4norm,A3),
  T is abs(A1*A2*A3),


  ((E = 0; T =< 0.5) -> Energy = 0 ;
   (E \= 0, T > 0.7) -> Energy = E*1000 ;
   (E \= 0, T > 0.5, T < 0.7) -> Energy = (E*1000*(T-0.5))/0.2 ).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%            LABELING HEURISTICS:         %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Tertiary = [ Ca list, Centroids list]

my_labeling(Code,Tertiary,ProbTot,Src,PrefDom,TotalEnergy,Torsion,ContactEnergy,DistanceList,ContactEnergyStructured):-
   retractall(best(_)), assert(best(0)),
   my_labeling_aux(Code,PrefDom,Tertiary,0,0,TotalEnergy,Torsion,ContactEnergy,DistanceList,ContactEnergyStructured),
   retractall(best(_)),
   assert(best(TotalEnergy)),

%%% dbg ellipsoid
%   ( %%%% OPTIONS - if then
%     options(ellipsoid(AACode, [RefX,RefY,RefZ], [[F1x,F1y,F1z],[F2x,F2y,F2z],Rad])),
%    write(test),
%    Tertiary=[CA,_],
%    ellipsoid_constraint(CA, AACode, [RefX,RefY,RefZ], [[F1x,F1y,F1z],[F2x,F2y,F2z],Rad]),
%    write(ok),nl
%     ;
%    \+options(ellipsoid(AACode, [RefX,RefY,RefZ], [[F1x,F1y,F1z],[F2x,F2y,F2z],Rad]))),



   final_data(Code,ProbTot,0,Src,1).

my_labeling_aux([],_,_,_,_,_,_,_,_,_).
my_labeling_aux([C|R],[Dom|PrefDom],Ter,N,SumPos,TotalEnergy,Torsion,ContactEnergy,DistanceList,ContactEnergyStructured):-

   (ground(C)   -> SumPosF=0;
    \+ground(C) -> SumPosF=1),
    member(C,Dom),  %%% choose C following decreasing probability ordering
   %indomain(C),   %%% choose C in ID increasing order
   sublist(Dom,[C],Pos,_,_),

   SumPos1 is SumPos + SumPosF * Pos,

   options(max_most_frequent(MMF)),
   options(max_non_optimal_choices(MNOC)),
   SumPosF * Pos <MMF, %% only MMF most frequent options
   SumPos1 <MNOC,      %% only MNOC choices non-most-frequent

   best(Best),
   (options(improve_energy), TotalEnergy #< Best;  %% bound the energy improving solutions
    \+options(improve_energy)
   ),
   cleantupla(_,_,_,_Freq,C,_ID),
   %% Debugging printing
   %write('### level '),
   %write(l(N,SumPos1)),write(' from '),write(ID), write(' '), write(code(C,'from',Dom)), nl,
   %fd_min(Torsion,FDmin),fd_max(Torsion,FDmax),fd_min(ContactEnergy,FDEmin),fd_max(ContactEnergy,FDEmax),fd_min(TotalEnergy,TotMin),fd_max(TotalEnergy,TotMax),
   %%write(b(total(TotMin,TotMax),tors(FDmin,FDmax),contact(FDEmin,FDEmax))),nl,
   %write_dom_list_mat(DistanceList),nl,
   %write(DistanceList),nl,
   %write_dom_list_mat_en(ContactEnergyStructured),nl,
   N1 is N+1,
   my_labeling_aux(R,PrefDom,Ter,N1,SumPos1,TotalEnergy,Torsion,ContactEnergy,DistanceList,ContactEnergyStructured).

final_data([],ProbTot,ProbSum,[],N) :-
   ProbTot is ProbSum/(N-1).
final_data([C|R],ProbTot,ProbSum,[PID|SRC],N) :-
   cleantupla(_,_,_,Freq,C,PID),  %%% The PID of the protein where the fragment is retrieved
   N1 is N+1,
   freq_factor(Fact),
   (Freq>0,!,
    ProbSum1 is ProbSum + log(10,(Freq/Fact));
    ProbSum1 is ProbSum
   ),
   final_data(R,ProbTot,ProbSum1,SRC,N1).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% LOCAL SEARCH MAIN AND MOVES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

choice([X|_], [P|_], Cumul, Rand, X) :-
    Rand < Cumul + P.
choice([_|Xs], [P|Ps], Cumul, Rand, Y) :-
    Cumul1 is Cumul + P,
    Rand >= Cumul1,
    choice(Xs, Ps, Cumul1, Rand, Y).
choice([X], [P], Cumul, Rand, X) :-
    Rand < Cumul + P.

choice(Xs, Ps, Y) :- random(R), choice(Xs, Ps, 0, R, Y).



local(ID,Code,Energy, Primary,Secondary,Tertiary,
             PD, TE,CE,DL,CEs,TCP,TCC):-
     \+ another_sol(ID,Code,Energy, Primary,Secondary,Tertiary,
             PD, TE,CE,DL,CEs,TCP,TCC),
     local(ID,Code,Energy, Primary,Secondary,Tertiary,
             PD, TE,CE,DL,CEs,TCP,TCC).

another_sol(ID,Code,Energy, _Primary, Secondary,Tertiary,
             PD, TE,CE,DL,CEs,TCP,TCC) :-
     last_sol(Lim,BestSol,BestTer),
      %%% Look for the first solution
     (BestSol = null ->
          write('Debug: Looking for the first solution'),nl,
          my_labeling(Code,Tertiary, _PCost,_Src,  PD, Energy, TE,CE,DL,CEs)
      ;
      BestSol \= null ->
          random(3,11,Type), %random(1,11,Type), %%% Type in 1..10
          write('Debug: random Type = '),write(Type),nl,
          (Type=<1 ->
             Lim1 is 5*Lim//6,
             Energy #> Lim1,
             length(BestSol,NB),
             random(1,NB,START),
             END is min(START+NB//2,NB), %%% Half variables free
             write('Debug: Trying Worsening sol, with free vars: '),
             write(START),write(' - '),write(END),nl,
             binds(1,START,END,Code,BestSol,Tertiary,BestTer),
             write('Debug: Looking for sol. with Energy >  '),write(Lim1),nl
          ;
          Type >2 ->
            Energy #< Lim,
            write('Debug: Looking for sol. with Energy < '), write(Lim),nl,

            %random(1,3,Move), %%% Move in 1..2

            %random(1,2,Move), %%% Only choose move = 1, and only do 1 type

            %choice([1,2],[0.8,0.2],Move), %%% 80 to large pivot

            choice([1,2],[0.2,0.8],Move), %%% 20 to large pivot

            (Move = 1 ->
               large_pivot(Code,BestSol,Tertiary,BestTer,Secondary)
            ;
             Move = 2 ->
              large_crankshaft(Code,BestSol,Tertiary,BestTer,Secondary)
            )
          ),
       Tertiary=[CA,_Centroids],
       append(Code,CA,CodeTertiary1),
       %(options(use_centroids),!,
       %append(CodeTertiary1,Centroids,CodeTertiary);
       CodeTertiary1=CodeTertiary
       ,
       %),

       term_variables(Code,_,VARLIBERE),
       length(VARLIBERE,KK),
       write('Debug: looking for a solution with '),
       write(KK), write(' free vars in Code'),nl,
%       time_out( labeling([], CodeTertiary), 60000, Flag),
       time_out( labeling([], CodeTertiary), 180000, Flag),
       (Flag == success  ->
            true
       ;
        Flag == time_out ->
            write('Warning: Time out in the labeling stage'),nl,
            fail
        )
     ),
     !,                                 %%% Il cut ferma la ricerca di
     retract(last_sol(_,_,_)),           %%% altre sol con lo stesso Lim
     retract(ct(N)),M is N + 1, assert(ct(M)),
     %%%%%%%% NEW: PATCH TO INCLUDE THE ORIENTATION ENERGY  February 10, 2011.
     Tertiary=[CAS,_],
     orientation_energy(CAS,ORE), OE is integer(ORE),
     TOT is Energy+OE,
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     format("Found solution No ~d with complete energy: ~d\n",[M,TOT]),

     %%% ASSERISCO IL VAORE SENZA LA COMPONENTE DI ORIENTAZIONE

     assert(last_sol(Energy,Code,Tertiary)),
     best([Val,_,_,_,_,_],_,_,_,_),
     (Val > TOT ->
         retract(best(_,Primary,_,_,_)),
         search_time(Time),
         assert(best([TOT,TE,CE,TCP,TCC,OE],Primary,Tertiary,Time,Code)),
         best_printing(ID);
      true),
     fail.

%%%% large_pivot: an interval is ND chosen.
%%% inner variables are "free", outer variables
%%% assume their previous value.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 08/12/2009

large_pivot(Vars,BestSol,Tertiary,BestTer,Pivots) :-
     length(Vars,M), N is M - 3,
     random_select(START,Pivots,_),
     A1 is min(N,START+3), % minimum interval = 4
     A3 is min(N,START+10),% maximum interval = 9
     random(A1,A3,END),
     format("Debug: PIVOT - free interval: ~d-~d \n",[START,END]),
     binds(1,START,END,Vars,BestSol,Tertiary,BestTer).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% scelgo due frammenti da 4 che lascio liberi
%%% Il primo pezzo di terziaria e' bloccato.
%%% In pratica si ruota dentro. L'inizio della
%%% parte finale e' fissato ma con tolleranza.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

large_crankshaft(Vars,BestSol,Tertiary,BestTer,Pivots)    :-
     length(Vars,M), N is M - 3,
     random_select(START1,Pivots,_),
     END1 is min(N,START1+3), % minimum interval = 4
     filter(END1,Pivots,REST),
     random_select(START2,REST,_),
     END2 is min(N,START2+3),
     format("Debug: CRANKSHAFT - free intervals ~d-~d and ~d-~d\n",[START1,END1, START2,END2]),
     binds(1,START1,END1,START2,END2,Vars,BestSol,Tertiary,BestTer).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% integer - constraint based - MATRIX ARITHMETIC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% identity_int returns a identity matrix with the
%%% allowed integer precision (e.g. F=1000)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%12/11/2009

identity_int([F,0,0,0,F,0,0,0,F]):-matrix_factor(F).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% matrix_mult_int computes the product ot two matrices %%%%%
%%% and rescale the result by matrix_factor F            %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%12/11/2009

matrix_mult_int([O11,O12,O13,O21,O22,O23,O31,O32,O33],
           [M11,M12,M13,M21,M22,M23,M31,M32,M33],
           [N11,N12,N13,N21,N22,N23,N31,N32,N33]) :-
     matrix_factor(F),
     N11 #= (O11 *M11 + O12 * M21 + O13 * M31)/F,
     N12 #= (O11 *M12 + O12 * M22 + O13 * M32)/F,
     N13 #= (O11 *M13 + O12 * M23 + O13 * M33)/F,
     N21 #= (O21 *M11 + O22 * M21 + O23 * M31)/F,
     N22 #= (O21 *M12 + O22 * M22 + O23 * M32)/F,
     N23 #= (O21 *M13 + O22 * M23 + O23 * M33)/F,
     N31 #= (O31 *M11 + O32 * M21 + O33 * M31)/F,
     N32 #= (O31 *M12 + O32 * M22 + O33 * M32)/F,
     N33 #= (O31 *M13 + O32 * M23 + O33 * M33)/F.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% vector_rot_int rotates a list of vectors Vx,Vy,Vz using the
%%%    rotation matrix [M11,...,M33] and
%%%    and rescale the result by matrix_factor F
%%% aux: vector_rot_int that rotates a single vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%12/11/2009

rotate_coord_int([],_,[]).
rotate_coord_int([X,Y,Z|R],MAT,[X1,Y1,Z1|S]) :-
    vector_rot_int([X,Y,Z],MAT,[X1,Y1,Z1]),
    rotate_coord_int( R,MAT,S).

vector_rot_int([V1,V2,V3],
            [M11,M12,M13, M21,M22,M23, M31,M32,M33],
            [N1,N2,N3]) :-
     matrix_factor(F),
     N1 #= (V1 *M11 + V2 * M21 + V3 * M31)/F,
     N2 #= (V1 *M12 + V2 * M22 + V3 * M32)/F,
     N3 #= (V1 *M13 + V2 * M23 + V3 * M33)/F.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% shift_int shifts a lists of points by DX,DY,DZ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%12/11/2009

shift_int([],_,_,_,[]).
shift_int([X,Y,Z|VIn],DX,DY,DZ,[Xo,Yo,Zo|VOut]) :-
    Xo #= X + DX,
    Yo #= Y + DY,
    Zo #= Z + DZ,
    shift_int(VIn,DX,DY,DZ,VOut).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% mult_all_int multiplies all points of a list by a Factor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%12/11/2009

mult_all_int([],_,[]).
mult_all_int([X|Rest], Factor,[X1|Rest1]) :-
    (ground(X),!,
    X1 is integer(X*Factor);
    true),
    mult_all_int(Rest, Factor, Rest1).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% floating point MATRIX ARITHMETIC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prod_vect executes the cross product of two 3D vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%12/11/2009

prod_vect(Ax,Ay,Az, Bx,By,Bz, Cx,Cy,Cz):-
   Cx is Ay*Bz-By*Az,
   Cy is Az*Bx-Bz*Ax,
   Cz is Ax*By-Bx*Ay.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prod_scal executes the scalar product of two 3D vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%12/11/2009

prod_scal(Ax,Ay,Az, Bx,By,Bz, S):-
  S is (Bx*Ax)+(By*Ay)+(Bz*Az).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% normalize scales a vector to a size "1"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%12/11/2009

normalize(X,Y,Z,Xn,Yn,Zn):-
   Norm is sqrt(X*X+Y*Y+Z*Z),
   (Norm>0,!, Xn is X/Norm, Yn is Y/Norm, Zn is Z/Norm;
    Norm=<0,!, Xn = 0, Yn = 0, Zn = 0).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%    AUXILIARY PREDICATES         %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% mutate/3 (aux mutate/4)
%%%% TO BE EXPLAINED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 08/02/2010

mutate(Primary,Primary,[]).
mutate(Primary,Primary1,[[Pos,NewAA]|List]):-
  mutate(Primary,Primary2,List),
  mutate(Primary2,Primary1,Pos,NewAA).

mutate(PrimaryIn,PrimaryOut,Pos,NewAA):-
    reference_aa_input_prot(Ofs),
    AANumber is Pos+1-Ofs,
    Pre is AANumber-1,
    sublist(PrimaryIn,PreList,0,Pre,_),
    append(PreList,[_Mut|Suf],PrimaryIn),
    append(PreList,[NewAA|Suf],PrimaryOut).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% randomize/0
%%%   sets random seeds using system date and time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

randomize:-
    datime(datime(A,B,C,D,E,F)),
    G is B*C,H is (D+1)*(E+1),
    F1 is F+1,
    setrand(random(F1,G,H,A)).

%%%%%%%% AUXILIARY predicates for local search:

filter(_A,[],[]).
filter(A,[B|R],R) :- A < B, !.
filter(A,[_|R],S) :-
     filter(A,R,S).

%%%%%%%%%%%%%%%%%%
%% patch/3 builds a list from two previous list.
%% If the element in the first is ground is selected,
%% otherwise the element in the second is chosen
%%%%%%%%%%%%%%%%%%

patch([],[],[]).
patch([A|R],[B|R1],[C|R2]):-
   (ground(A),!, C=A ;
    C=B ),
   patch(R,R1,R2).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% aux:
%%%  binds/7 retains the bindings in the not free
%%%  parts and in th efirst part of the tertiary.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

binds(C,A,B,[V|Ars],[V|Est],[[X,Y,Z|T],[XC,YC,ZC|TC]],[[X,Y,Z|R],[XC,YC,ZC|RC]]) :-
     C < A,  C1 is C + 1,
     binds(C1,A,B,Ars,Est,[T,TC],[R,RC]).

binds(C,A,B,[_V|Ars],[_B|Est],_,_) :-
     C >= A, C =< B, C1 is C + 1,
     binds(C1,A,B,Ars,Est,_,_).

binds(C,_A,B,Best,Best,_,_) :-
     C > B.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% aux: binds/9 sets the variables for the crashaft move
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

binds(C,A,A1,B,B1,[V|Ars],[V|Est],[[X,Y,Z|T],[XC,YC,ZC|TC]],[[X,Y,Z|R],[XC,YC,ZC|RC]]) :-
     C < A,  C1 is C + 1,
     binds(C1,A,A1,B,B1,Ars,Est,[T,TC],[R,RC]).

binds(C,A,A1,B,B1,[_V|Ars],[_B|Est],[[_,_,_|T],[_,_,_|TC]],[[_,_,_|R],[_,_,_|RC]]) :-
     C >= A, C =< A1, C1 is C + 1,
     binds(C1,A,A1,B,B1,Ars,Est,[T,TC],[R,RC]).

binds(C,A,A1,B,B1,[V|Ars],[V|Est],[[_,_,_|T],[_,_,_|TC]],[[_,_,_|R],[_,_,_|RC]]) :-
     C > A1, C < B, C1 is C + 1,
     binds(C1,A,A1,B,B1,Ars,Est,[T,TC],[R,RC]).

binds(C,A,A1,B,B1,[_V|Ars],[_B|Est],[[_,_,_|T],[_,_,_|TC]],[[_,_,_|R],[_,_,_|RC]]) :-
     C >= B,  C =< B1, C1 is C + 1,
     binds(C1,A,A1,B,B1,Ars,Est,[T,TC],[R,RC]).

binds(C,_A,_A1,_B,B1,Best,Best,[[X1,Y1,Z1|_],[_,_,_|_]],[[XT,YT,ZT|_],[_,_,_|_]]) :-
     C > B1,
     BOUND = 1000,
     abs(X1 - XT) #< BOUND,
     abs(Y1 - YT) #< BOUND,
     abs(Z1 - ZT) #< BOUND.


%%% aux: it selects points not involved in
%%% secondary predictions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

good_points(M,N,_,[]) :- M > N,!.
good_points(M,N,Secondary,R) :-
     (member(helix(A,B),Secondary);
      member(strand(A,B),Secondary)),
      M >= A,
      M =< B,
      !,
      M1 is M + 1,
      good_points(M1,N,Secondary,R).
good_points(M,N,Secondary,[M|R]) :-
      M1 is M + 1,
      good_points(M1,N,Secondary,R).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% diameter computes the known diameter of a predicted
%%%% aux: diameter_aux computes the max distance of a point to the others
%%%% the diameter is computed on CA distances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 13/11/09

diameter([CA,_],D):-
  diameter1(CA,D).

diameter1([_,_,_],0).
diameter1([X,Y,Z|R],Max):-
   diameter_aux(R,X,Y,Z,Max1),
   diameter1(R,Max2),
   (Max1>Max2,!,Max=Max1;
    Max=Max2).

diameter_aux([],_X,_Y,_Z,0).
diameter_aux([X1,Y1,Z1|R],X,Y,Z,Max):-
   resolution(Res),
   Dist is sqrt(((X-X1)*(X-X1) + (Y-Y1)*(Y-Y1) + (Z-Z1)*(Z-Z1))/(Res*Res)),
   diameter_aux(R,X,Y,Z,Max1),
   (Dist > Max1,!,Max=Dist;
    Max=Max1).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% reset cleans the workspace and stores the initial runtime
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 13/11/2009

reset :-
   retractall(s(_)),
   retractall(ct(_)), assert(ct(0)),
   retractall(time(_,_)),
   statistics(runtime,[T,_]),
   assert(time(initial,T)).

new_options(Opt):-
   retractall(options(_)),
   check_options(Opt),
   set_options(Opt),

   %% if not specified, add default (+inf)
   (options(max_most_frequent(_)),!;
    assert(options(max_most_frequent(10000)))),

   (options(max_non_optimal_choices(_)),!;
   assert(options(max_non_optimal_choices(10000)))),

  write('The following options are being used: \n'),
  write_options1.

set_options([]).
set_options([O|R]):-
  assert(options(O)), set_options(R).


write_options:-
  write('Change them by calling the predicate: \n\tnew_options([opt1, opt2, ...]).\n\n'),
  write('Available options:\n'),
  write_option.

write_option:-
  available_option(Opt,Desc),
  write('\t '),write(Opt),write(': '),write(Desc),nl,
  fail.
write_option.


check_options([]).
check_options([O|R]):-
   (available_option(O,_),!;
    write('Error: option '),write(O),write(' not available'),nl,!,
    fail
   ),
   check_options(R).

available_option(use_secondary ,'impose the secondary structure contained in prot_list.txt file').
available_option(use_centroids ,'use a model with centroids (otherwise only CA backbone)').
available_option(improve_energy , 'every time a structure is found, look for a better energy (otherwise simple conformational search)').
available_option(use_original_codes ,'overimpose the original templates everywhere').
available_option(use_original_codes(_Ranges) ,'use the original torsions from templates on the specified Ranges').
available_option(use_original_tertiary(_Ranges) ,'fix the original tertiary position for the specified Ranges').
available_option(mutation(_Mutation) ,'mutates the input sequence with the list of pairs [Position,NewAminoAcid1Letter]').
available_option(box(_AA, [_X,_Y,_Z], [[_,_],[_,_],[_,_]]) ,'box domain: AA number, X,Y,Z original ref for first AA, minmax x y z').
available_option(ellipsoid(_AA, [_X,_Y,_Z], [[_,_,_],[_,_,_],_]) ,'box domain: AA number, X,Y,Z original ref for first AA, F1 [x,y,z], F2 [x,y,z], radius').
available_option(contiguous_check ,'perform a check for ca ca distance of 3.6--4.0 for each structure').
available_option(max_most_frequent(_MMF) ,'MMF = number of most frequent fragments').
available_option(max_non_optimal_choices(_MNOC) ,'MNOC = number of non total non-most-frequent choices').

write_options1:-
   options(A),
   write('\t'),write(A),nl,
   fail.
write_options1:-nl.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% loading_time writes the loading time and stores its
%%% absolute value.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 13/11/2009

loading_time :-
   statistics(runtime,[Tl,_]),
   assert(time(loading,Tl)),
   time(initial,Ti),
   Time is (Tl - Ti)/1000,
   format("### Data processed and loaded in ~3F s\n",[Time]).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% constraint_time writes the loading time and stores its
%%% absolute value.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 13/11/2009

constraint_time  :-
   statistics(runtime,[Tl,_]),
   assert(time(constraint,Tl)),
   time(loading,Ti),
   Time is (Tl - Ti)/1000,
   format("### Constraints added in ~3F s\n",[Time]).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% search_time computes the (total) search time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 13/11/2009

search_time(Time) :-
   statistics(runtime,[Tl,_]),
   time(constraint,Ti),
   Time is (Tl - Ti)/1000.
   %, format("### Search time:  ~3F s\n",[Time]).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% proteina retrieves the information of a known
%%% or of an unknown protein
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%16/11/2009

proteina(ID,Primary,Secondary,Reference_aa_input_prot) :-
    protein(ID,Primary,Secondary,Reference_aa_input_prot),!,
    retractall(reference_aa_input_prot(_)),
    assert(reference_aa_input_prot(Reference_aa_input_prot)).
proteina(ID,Primary,Secondary,Reference_aa_input_prot) :-
    new_protein(ID,Primary,Secondary,Reference_aa_input_prot),
    retractall(reference_aa_input_prot(_)),
    assert(reference_aa_input_prot(Reference_aa_input_prot)).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% split_matrix/3 extractw two flat lists (non contiguous
%%% and contiguous) from the triangular distances matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%14/11/2009

split_matrix([],[],[]).
split_matrix([[_,[]]|_],[],[]). %% maybe useless
split_matrix([[_,[D|Ds]]|R],FlatALL,[D|NEXTs]) :-
     split_matrix(R , FlatREST,NEXTs),
     append(Ds,FlatREST,FlatALL).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% rev_sort order the list in descending order
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%16/11/2009

rev_sort(List, RevSortedList) :-
      sort(List,SortedList),
      reverse(SortedList,RevSortedList).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  same_value(LIST,VAL) assign VAL to all
%%%  elements of LIST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 14/11/2009

same_value([],_).
same_value([A|R],A):-
  same_value(R,A).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  same_constraint(LIST,OP,VAL) forces the constraint
%%%  X op VAL to all elements in LIST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 14/11/2009

same_constraint([],_,_).
same_constraint([X|LIST],OP,VAL):-
  C =.. [OP,X,VAL],
  C,
  same_constraint(LIST,OP,VAL).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  contiguous_check check
%%%    prints a warning if, due to approximation errors,
%%%    the roughly 3.8AA are not satisfied
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%13/11/2009

contiguous_check([_,_,_],_).
contiguous_check([X1,Y1,Z1,X2,Y2,Z2|Tertiary],N):-
   N1 is N+1,
   resolution(Res),
   Delta is sqrt((X2-X1)*(X2-X1)+(Y2-Y1)*(Y2-Y1)+(Z2-Z1)*(Z2-Z1))/Res,
   Err is abs(Delta - 3.81),
  (Err>0.2 -> format("!!! Warning. Consecutive aminoacids ~d-~d+1 are far ~3F AA\n",[N,N1,Delta])
   ;
   Err=<0.2 -> true),
   contiguous_check([X2,Y2,Z2|Tertiary],N1).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

convert_list([],[]).
convert_list([A|R],[CA|CR]):-
    name_conversion(A3, A, _),
    class(A3,CA),
    convert_list(R,CR).
convert_list_temp([],[]).
convert_list_temp([A|R],[CA|CR]):-
    name_conversion(A3, A, _),
    temp(class(A3,CA)),
    convert_list_temp(R,CR).

%%% convert the template positions to FD according to the
%%% resolution
cleantupla(Amino,VFD,VFDC,Prob,Id,Prot) :-
   tuple(Amino,V,[C2,C3],Prob,Id,Prot),
   resolution(R),
   mult_all_int(V,R,VFD),
   rec_centr_scale(C2,VC2,R),
   rec_centr_scale(C3,VC3,R),
   VFDC=[VC2,VC3].
rec_centr_scale([],[],_).
rec_centr_scale([[AA,Point]|R],[[AA,Point1]|R1],Sc):-
   mult_all_int(Point,Sc,Point1),
   rec_centr_scale(R,R1,Sc).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% STYLE FOR INPUT DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% A set of four/five-tuplets
%%% tupla(Potlist,[PosList,CPosList],Prob,Id,OrigProt)
%%% where Potlis is the list of amino acid names,
%%% PosList is the positions of the Calpha
%%% PosList is the positions of the Centroids (cbeta)
%%% Prob is the probability of that local conformation
%%% Id is the tuple identifier.
%%% OrigProt is the original protein that provides this template.

%%% The predicate   next(ID1,ID2,M)
%%%%%% holds if the sequence 2-5 of ID1 and 1-4 are the same
%%%%%% tupla([_,B,C,D,E],
%%%%%% [_,_,_,X2,Y2,Z2,...,X5,Y5,Z5],_,ID1,PID),
%%%%%% tupla([B,C,D,E,_],
%%%%%% [Xa1,Ya1,Za1,...,Xa4,Ya4,Za4,_,_,_],_,ID2,PID),
%%%%%% Moreover, the RMSD of
%%%%%% [X2,Y2,Z2,...,X5,Y5,Z5] and
%%%%%% [Xa1,Ya1,Za1,...,Xa4,Ya4,Za4]
%%%%%% is below a give threshold (to be chosen)
%%%%%% M is the rotation matrix to align the second on the first.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% FILES CONSULTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

:-compile('tables.pl').
:-compile('prot-list.pl').
:-compile('prot_preproc.pl').
:-compile('printing_predicates.pl').

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Resolution (how many times an Angstrom is divided)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

resolution(100).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Coefficient for rotation matrices (they are stored as integers)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matrix_factor(1000).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Coefficient for frequency (they are stored as integers)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq_factor(1000).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Length of input data: currently 4-tuplets or 5-tuplets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datalength(4).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Search parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

min_distance(3.2).     %%% between ca -ca and ca -cg
min_distance_cg(1.0).  %%% between cg -cg
max_distance(100).     %%%
distanceSsbond(6.6).   %%% max distance between two Ca connected by an ssbond

toleranceLinkDist(1.5).  %%% when connecting a CA to a already fixed substructure (in space), allow X Angstrom in the link distance

%%% used for estimation of max distance between aa i and j (= abs(i-j)/2*distanceCaCaStep2+abs(i-j)%2*distanceCaCa)
distanceCaCa(3.9).       %%% max distance between two consecutive Ca (if decreased could cut out some templates -> the search fail!)
distanceCaCaStep2(7.5).  %%% max distance between two Ca at distance 2 (limited by the widest bend angle) (if decreased could cut out some templates -> the search fail!)
next_CaCadistance(3.8).  %% exact distance between consecutive aa (imposed as a constrain with a +-5% of error)

:-write('Set default options\n'),
  new_options([use_secondary,improve_energy,contiguous_check,use_centroids]),
%%  new_options([use_original_codes,use_centroids]),
%%  new_options([use_secondary,contiguous_check,use_centroids]),
  write_options.

:- randomize.

%new_options([use_secondary,use_original_codes([[307,523]]),box(540, [9.931,  12.549,   9.590], [[25,50],[-28,2],[-4,30]])]).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%     END of FILE    %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
