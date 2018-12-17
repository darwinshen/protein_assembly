%%% Some TESTED DATA
%%% All proteins and their secondary structure information are
%%% picked from the Protein Data Bank. Their ID is the same used
%%% by PDB. 
%%% Note: first amino acid is at position 1 in the list if protein/3 is used
%%% Last modification: 091218
 
protein(ID, Primary, Secondary, 1):-
  protein(ID, Primary, Secondary).
  
  
protein(ID, Primary, Secondary, FirstAaPosition):-
     ID = '1XPC',
	 FirstAaPosition = 307,
     Primary = 
[a,l,s,l,t,a,d,q,m,v,
s,a,l,l,d,a,e,p,p,i,
l,y,s,e,y,d,p,t,r,p,
f,s,e,a,s,m,m,g,l,l,
t,n,l,a,d,r,e,l,v,h,
m,i,n,w,a,k,r,v,p,g,
f,v,d,l,t,l,h,d,q,v,
h,l,l,e,c,a,w,l,e,i,
l,m,i,g,l,v,w,r,s,m,
e,h,p,g,k,l,l,f,a,p,
n,l,l,l,d,r,n,q,g,k,
c,v,e,g,m,v,e,i,f,d,
m,l,l,a,t,s,s,r,f,r,
m,m,n,l,q,g,e,e,f,v,
c,l,k,s,i,i,l,l,n,s,
g,v,y,t,f,l,s,s,t,l,
k,s,l,e,e,k,d,h,i,h,
r,v,l,d,k,i,t,d,t,l,
i,h,l,m,a,k,a,g,l,t,
l,q,q,q,h,q,r,l,a,q,
l,l,l,i,l,s,h,i,r,h,
m,s,n,k,g,m,e,h,l,y,
s,m,k,c,k,n,v,v,p,l,
y,d,l,l,l,e,m,l,d,a
%,h,r,l,h,a
],
     Secondary = [ helix(538,546) ].

protein(ID, Primary, Secondary):-
     ID = 'TEST',
     Primary = [c,c,c,c] ,
     Secondary = [  ].

protein(ID, Primary, Secondary):-
%% length: 12
     ID = '1LE0',
     Primary = [s,w,t,w,e,g,n,k,w,t,w,k] ,
     Secondary = [ strand(2,5), strand(8,11) ].


protein(ID, Primary, Secondary):-
%% length: 12
        ID = '1KVG',
        Primary = [s,c,h,f,g,p,l,g,w,v,c,k],
        Secondary = [ ssbond(2,11),
                      strand(3,4), strand(9,11)].

protein(ID, Primary, Secondary):-
%% length: 12
        ID = '1KVGnosec',
        Primary = [s,c,h,f,g,p,l,g,w,v,c,k],
        Secondary = [sicstus].


protein(ID, Primary, Secondary):-
%%% Note: cut an "x" in the last position.
%% length: 16
         ID = '1LE3',   
       Primary = [g,e,w,t,w,d,d,a,t,k,
                    t,w,t,w,t,e], 
       Secondary = [strand(2,6), strand(11,15)]. 


protein(ID, Primary, Secondary):-
%%% Note: cut an "x" in the last position.
%% length: 18
        ID = '1PG1',
       Primary = [r,g,g,r,l,c,y,c,r,r,
                   r,f,c,v,c,v,g,r], 
      Secondary = [ssbond(6,15), ssbond(8,13),
                     strand(4,9), strand(12,17) ].
        
protein(ID, Primary, Secondary):-
%%% length: 17
        ID = '1EDP',
         Primary = [c,s,c,s,s,l,m,d,k,e,
                    c,v,y,f,c,h,l],
         Secondary = [ssbond(1,15), ssbond(3,11),
                      helix(9,15) 
                     ].
        
protein(ID, Primary, Secondary):-
%%% Note: cut an "x" in the last position.
%%% length: 34  ---
        ID = '1ZDD',
       Primary = [f,n,m,q,c,q,r,r,f,y,
                  e,a,l,h,d,p,n,l,n,e,
                  e,q,r,n,a,k,i,k,s,i,
                  r,d,d,c],
        Secondary = [
                 ssbond(5, 34),
                      helix(3,13),helix(20,33)].


protein(ID, Primary, Secondary):-
%%% length: 36
       ID = '1VII',
        Primary = [m,l,s,d,e,d,f,k,a,v, 
                   f,g,m,t,r,s,a,f,a,n, 
                      l,p,l,w,k,q,q,n,l,k, 
                      k,e,k,g,l,f],
       Secondary = [helix(4,8),helix(15,18),helix(23,32)].


protein(ID, Primary, Secondary):-
%%% length: 36
       ID = '1VIIgoriv',
        Primary = [m,l,s,d,e,d,f,k,a,v, 
                   f,g,m,t,r,s,a,f,a,n, 
                      l,p,l,w,k,q,q,n,l,k, 
                      k,e,k,g,l,f],
       Secondary = [ %%% BY GORIV
       helix(8,17),
       helix(21,30),
       strand(34,35)].


protein(ID, Primary, Secondary):-
%%% length: 37
        ID = '1E0M',
         Primary = [s,m,g,l,p,p,g,w,d,e, 
                    y,k,t,h,n,g,k,t,y,y, 
                        y,n,h,n,t,k,t,s,t,w, 
                        t,d,p,r,m,s,s],
         Secondary = [strand(8,12),strand(18,22),strand(27,29)].


        
protein(ID, Primary, Secondary):-
%%% length: 40
        ID = '2GP8',
       Primary = [i,t,g,d,v,s,a,a,n,k, 
                    d,a,i,r,k,q,m,d,a,a, 
                        a,s,k,g,d,v,e,t,y,r, 
                        k,l,k,a,k,l,k,g,i,r],
        Secondary = [helix(6,21),helix(26,38)].     


protein(ID, Primary, Secondary):-
%% length: 46
     ID = '1ED0',
     Primary = [k,s,c,c,p,n,t,t,g,r,
                n,i,y,n,a,c,r,l,t,g,
                a,p,r,p,t,c,a,k,l,s,
                g,c,k,i,i,s,g,s,t,c,
                p,s,d,y,p,k],
     Secondary = [
                   helix(6,18),helix(22,30),
                   strand(1,5),strand(31,34)].


protein(ID, Primary, Secondary, FirstAaPosition):-
%% length: 54
     ID = '1ENH',
	 FirstAaPosition = 3,
     Primary = [r,p,r,t,a,f,s,s,e,q,
                 l,a,r,l,k,r,e,f,n,e,
                     n,r,y,l,t,e,r,r,r,q,
                     q,l,s,s,e,l,g,l,n,e,
                     a,q,i,k,i,w,f,q,n,k,
                     r,a,k,i],
     Secondary = [helix(10,22),helix(28,38),helix(42,56)].
protein(ID, Primary, Secondary):-
%% length: 54
     ID = '1ENHweak',
     Primary = [ r,p,r,t,a,f,s,s,e,q,
                 l,a,r,l,k,r,e,f,n,e,
                 n,r,y,l,t,e,r,r,r,q,
                 q,l,s,s,e,l,g,l,n,e,
                 a,q,i,k,i,w,f,q,n,k,
                 r,a,k,i ],
     Secondary = [helix(8,18),helix(28,36),helix(40,52)].
       
protein(ID, Primary, Secondary):-
%%% length = 58
      ID = '6PTI',
      Primary = [r,p,d,f,c,l,e,p,p,y, 
                 t,g,p,c,k,a,r,i,i,r, 
                     y,f,y,n,a,k,a,g,l,c, 
                     q,t,f,v,y,g,g,c,r,a, 
                     k,r,n,n,f,k,s,a,e,d, 
                     c,m,r,t,c,g,g,a],
      Secondary = [ssbond(5,55),ssbond(14,38),
                        ssbond(30,51),
                        helix(3,6), 
                        helix(48,55),
                        strand(18,24), 
                        strand(29,35)].


protein(ID, Primary, Secondary):-
%%% length = 60
  ID = '2IGD',
  Primary = [m,t,p,a,v,t,t,y,k,l,
             v,i,n,g,k,t,l,k,g,e,
             t,t,t,k,a,v,d,a,e,t,
             a,e,k,a,f,k,q,y,a,n,
             d,n,g,v,d,g,v,w,t,y,
             d,d,a,t,k,t,f,t,v,t],
  Secondary = [ helix(28,41), strand(18,25), strand(6,13), 
                strand(56,60), 
                strand(47,51)].


protein(ID, Primary, Secondary):-
%%% length = 61
  ID = '2ERA',
  Primary = [r,i,c,f,n,h,q,g,s,q,
             p,q,t,t,k,t,c,s,p,g,
             e,s,s,c,y,n,k,q,w,s,
             d,f,r,g,t,i,i,e,r,g,
             c,g,c,p,t,v,k,p,g,i,
             k,l,s,c,c,e,s,e,v,c,
             n],
  Secondary = [ strand(2,5), strand(13,17), strand(34,41), 
                strand(23,31), strand(50,55)
			  ].



protein(ID, Primary, Secondary):-
%%% length = 63
  ID = '1SN1',
  Primary = [v,r,d,a,y,i,a,k,p,h,
             n,c,v,y,e,c,a,r,n,e,
             y,c,n,d,l,c,t,k,n,g,
             a,k,s,g,y,c,q,w,v,g,
             k,y,g,n,g,c,w,c,i,e,
             l,p,d,n,v,p,i,r,v,p,
             g,k,c],
  Secondary = [ helix(18,28), strand(2,6), strand(32,39), strand(43,51)].


protein(ID, Primary, Secondary):-
%%% length = 63
  ID = '1YPA',
  Primary = [m,k,t,e,w,p,e,l,v,g,
             k,a,v,a,a,a,k,k,v,i,
             l,q,d,k,p,e,a,q,i,i,
             v,l,p,v,g,t,i,v,t,m,
             e,y,r,i,d,r,v,r,l,f,
             v,d,k,l,d,n,i,a,q,v,
             p,r,v],
  Secondary = [ helix(31,43), strand(24,29)].
       
	   
protein(ID, Primary, Secondary, FirstAaPosition):-
%% length: 54
     ID = '2KAP',
	 FirstAaPosition = 1,
     Primary = [k,e,a,c,d,w,l,r,a,t,
g,f,p,q,y,a,q,l,y,e,
d,f,l,f,p,i,d,i,s,l,
v,k,r,e,h,d,f,l,d,r,
d,a,i,e,a,l,c,r,r,l,
n,t,l,n,k,c,a,v,m,k],
     Secondary = [helix(1,10),
	              helix(14,21),helix(27,35),helix(40,58)].

protein(ID, Primary, Secondary, FirstAaPosition):-
%% length: 54
     ID = '2K9D',
	 FirstAaPosition = 462,
     Primary = [
s,v,i,r,s,i,i,k,s,s,
r,l,e,e,d,r,k,r,y,l,
m,t,l,l,d,d,i,k,g,a,
n,d,l,a,k,f,h,q,m,l,
v,k,i,i],
     Secondary = [helix(462,471),helix(474,486),helix(490,505)].

	 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

protein(ID, Primary, Secondary):-
%% protein length: 69
  ID = '1AIL',
  Primary = [m,d,s,n,t,v,s,s,f,q,
             v,d,c,f,l,w,h,v,r,k,
             q,v,v,d,q,e,l,g,d,a,
             p,f,l,d,r,l,r,r,d,q,
             k,s,l,r,g,r,g,s,t,l,
             g,l,n,i,e,a,a,t,h,v,
             g,k,q,i,v,e,k,i,l],
  Secondary = [ helix(3,25), helix(30,50), helix(54,68)].


protein(ID, Primary, Secondary):-
%% protein length: 71
  ID = '1FVS',
  Primary = [a,r,e,v,i,l,a,v,h,g,
             m,t,c,s,a,c,t,n,t,i,
             n,t,q,l,r,a,l,k,g,v,
             t,k,c,d,i,s,l,v,t,n,
             e,c,q,v,t,y,d,n,e,v,
             t,a,d,s,i,k,e,i,i,e,
             d,c,g,f,d,c,e,i,l,r,
             d],
  Secondary = [ helix(13,26), helix(51,62), strand(30,36), strand(41,46), strand(2,9), strand(65,70)].


protein(ID, Primary, Secondary):-
%% protein length: 71
  ID = '1AW0',
  Primary = [l,t,q,e,t,v,i,n,i,d,
             g,m,t,c,n,s,c,v,q,s,
             i,e,g,v,i,s,k,k,p,g,
             v,k,s,i,r,v,s,l,a,n,
             s,n,g,t,v,e,y,d,p,l,
             l,t,s,p,e,t,l,r,g,a,
             i,e,d,m,g,f,d,a,t,l,
             s],
  Secondary = [ helix(16,26), helix(54,64), strand(31,37), strand(42,47), strand(3,10), strand(67,70)].


protein(ID, Primary, Secondary):-
%% protein length: 72
  ID = '1FD8',
  Primary = [m,a,e,i,k,h,y,q,f,n,
             v,v,m,t,c,s,g,c,s,g,
             a,v,n,k,v,l,t,k,l,e,
             p,d,v,s,k,i,d,i,s,l,
             e,k,q,l,v,d,v,y,t,t,
             l,p,y,d,f,i,l,e,k,i,
             k,k,t,g,k,e,v,r,s,g,
             k,q],
  Secondary = [ helix(17,30), helix(52,63), strand(35,39), strand(44,49), strand(5,11), strand(67,73)].




protein(ID, Primary, Secondary):-
%% protein length: 78
  ID = '1L6T',
  Primary = [m,e,n,l,n,m,d,l,l,y,
             m,a,a,a,v,m,m,g,l,a,
             a,i,g,d,a,i,g,i,g,i,
             l,g,g,k,f,l,e,g,a,a,
             r,q,p,d,l,i,p,l,l,r,
             t,q,f,f,i,v,m,g,l,v,
             n,a,i,p,m,i,a,v,g,l,
             g,l,y,v,m,f,a,v],
  Secondary = [ helix(4,40), helix(45,77)].



protein(ID, Primary, Secondary):-
%% protein length: 79
  ID = '1NH9',
  Primary = [m,d,n,v,v,l,i,g,k,k,
             p,v,m,n,y,v,v,a,v,l,
             t,q,l,t,s,n,d,e,v,i,
             i,k,a,r,g,k,a,i,n,k,
             a,v,d,v,a,e,m,i,r,n,
             r,f,i,k,d,i,k,i,k,k,
             i,e,i,g,t,d,k,e,v,n,
             v,s,t,i,e,i,v,l,a],
  Secondary = [ helix(11,26), helix(37,52), strand(4,7), strand(28,34), strand(78,86), strand(57,66)].


protein(ID, Primary, Secondary):-
%% protein length: 79
  ID = '1B4R',
  Primary = [a,t,l,v,g,p,h,g,p,l,
             a,s,g,q,l,a,a,f,h,i,
             a,a,p,l,p,v,t,a,t,r,
             w,d,f,g,d,g,s,a,e,v,
             d,a,a,g,p,a,a,s,h,r,
             y,v,l,p,g,r,y,h,v,t,
             a,v,l,a,l,g,a,g,s,a,
             l,l,g,t,d,v,q,v,e],
  Secondary = [ strand(2,4), strand(14,21), strand(46,51), strand(39,43), strand(28,32), strand(55,64), strand(69,78)].


protein(ID, Primary, Secondary):-
%% protein length: 83
  ID = '1W53',
  Primary = [m,d,f,r,e,v,i,e,q,r,
             y,h,q,l,l,s,r,y,i,a,
             e,l,t,e,t,s,l,y,q,a,
             q,k,f,s,r,k,t,i,e,h,
             q,i,p,p,e,e,i,i,s,i,
             h,r,k,v,l,k,e,l,y,p,
             s,l,p,e,d,v,f,h,s,l,
             d,f,l,i,e,v,m,i,g,y,
             g,m,a],
  Secondary = [ helix(2,22), helix(23,40), helix(43,59), helix(63,84)].


protein(ID, Primary, Secondary):-
%% protein length: 83
  ID = '1K1C',
  Primary = [v,q,q,k,v,e,v,r,l,k,
             t,g,l,q,a,r,p,a,a,l,
             f,v,q,e,a,n,r,f,t,s,
             d,v,f,l,e,k,d,g,k,k,
             v,n,a,k,s,i,m,g,l,m,
             s,l,a,v,s,t,g,t,e,v,
             t,l,i,a,q,g,e,d,e,q,
             e,a,l,e,k,l,a,a,y,v,
             q,e,e],
  Secondary = [ helix(15,28), helix(45,50), helix(68,82), strand(2,8), strand(59,66), strand(30,36), strand(39,41)].


protein(ID, Primary, Secondary):-
%% protein length: 86
  ID = '1NFJ',
  Primary = [e,h,v,v,y,v,g,n,k,p,
             v,m,n,y,v,l,a,t,l,t,
             q,l,n,e,g,a,d,e,v,v,
             i,k,a,r,g,r,a,i,s,r,
             a,v,d,v,a,e,i,v,r,n,
             r,f,m,p,g,v,k,v,k,e,
             i,k,i,d,t,e,e,l,e,s,
             e,q,g,r,r,s,n,v,s,t,
             i,e,i,v,l,a],
  Secondary = [ helix(10,24), helix(36,52), strand(2,5), strand(28,34), strand(77,86), strand(57,67)].


protein(ID, Primary, Secondary):-
%% protein length: 87
  ID = '1TIG',
  Primary = [i,n,v,k,e,v,r,l,s,p,
             t,i,e,e,h,d,f,n,t,k,
             l,r,n,a,r,k,f,l,e,k,
             g,d,k,v,k,a,t,i,r,f,
             k,g,r,a,i,t,h,k,e,i,
             g,q,r,v,l,d,r,l,s,e,
             a,c,a,d,i,a,v,v,e,t,
             a,p,k,m,d,g,r,n,m,f,
             l,v,l,a,p,k,n],
  Secondary = [ helix(14,29), helix(47,61), strand(3,8), strand(33,39), strand(78,85), strand(66,70)].



protein(ID, Primary, Secondary):-
%% protein length: 91
  ID = '3NCM',
  Primary = [y,v,m,f,k,n,a,p,t,p,
             q,e,f,k,e,g,e,d,a,v,
             i,v,c,d,v,v,s,s,l,p,
             p,t,i,i,w,k,h,k,g,r,
             d,v,i,l,k,k,d,v,r,f,
             i,v,l,s,n,n,y,l,q,i,
             r,g,i,k,k,t,d,e,g,t,
             y,r,c,e,g,r,i,l,a,r,
             g,e,i,n,f,k,d,i,q,v,
             i],
  Secondary = [ ssbond(23,73), helix(42,45), helix(64,67), strand(2,5), strand(24,27), strand(11,12), strand(82,91), strand(69,77), strand(32,36), strand(17,21), strand(58,62), strand(50,52)].




protein(ID, Primary, Secondary):-
%% protein length: 91
  ID = '1JNS',
  Primary = [a,k,t,a,a,a,l,h,i,l,
             v,k,e,e,k,l,a,l,d,l,
             l,e,q,i,k,n,g,a,d,f,
             g,k,l,a,k,k,h,s,i,c,
             p,s,g,k,r,g,g,d,l,g,
             e,f,r,q,g,q,m,v,p,a,
             f,d,k,v,v,f,s,c,p,v,
             l,e,p,t,g,p,l,h,t,q,
             f,g,y,h,i,i,k,v,l,y,
             r],
  Secondary = [ helix(14,27), helix(30,36), helix(43,47), helix(60,67), strand(50,53), strand(3,11), strand(82,89), strand(74,79)].


protein(ID, Primary, Secondary):-
%% protein length: 92
  ID = '1WIT',
  Primary = [l,k,p,k,i,l,t,a,s,r,
             k,i,k,i,k,a,g,f,t,h,
             n,l,e,v,d,f,i,g,a,p,
             d,p,t,a,t,w,t,v,g,d,
             s,g,a,a,l,a,p,e,l,l,
             v,d,a,k,s,s,t,t,s,i,
             f,f,p,s,a,k,r,a,d,s,
             g,n,y,k,l,k,v,k,n,e,
             l,g,e,d,e,a,i,f,e,v,
             i,v],
  Secondary = [ helix(67,69), strand(4,6), strand(20,27), strand(57,62), strand(50,54), strand(10,14), strand(82,93), strand(71,79), strand(33,37), strand(44,45)].



protein(ID, Primary, Secondary):-
%% protein length: 94
  ID = '1UFW',
  Primary = [g,s,s,g,s,s,g,s,s,f,
             q,g,p,l,d,a,t,v,v,v,
             n,l,q,s,p,t,l,e,e,k,
             n,e,f,p,e,d,l,r,t,e,
             l,m,q,t,l,g,s,y,g,t,
             i,v,l,v,r,i,n,q,g,q,
             m,l,v,t,f,a,d,s,h,s,
             a,l,s,v,l,d,v,d,g,m,
             k,v,k,g,r,a,v,k,i,s,
             g,p,s,s],
  Secondary = [ helix(4,7), helix(27,31), helix(35,48), helix(69,77), strand(53,57), strand(60,64), strand(17,22), strand(85,90), strand(81,82)].



protein(ID, Primary, Secondary):-
%% protein length: 94
  ID = '1CI5',
  Primary = [s,s,q,q,i,y,g,v,k,y,
             g,n,v,t,f,h,v,p,s,n,
             q,p,l,k,e,v,l,w,k,k,
             q,k,d,k,v,a,e,l,e,n,
             s,e,f,r,a,f,s,s,f,k,
             n,r,v,y,l,d,t,k,s,g,
             s,l,t,i,y,n,l,t,s,s,
             d,e,d,e,y,e,m,e,s,p,
             n,i,t,d,s,m,k,f,f,l,
             y,v,g,e],
  Secondary = [ helix(49,52), strand(3,8), strand(85,93), 
                strand(74,79), strand(26,30), strand(33,39), strand(42,45)].


protein(ID, Primary, Secondary):-
%% protein length: 95
  ID = '1LWR',
  Primary = [a,g,p,s,a,p,k,l,e,g,
             q,m,g,e,d,g,n,s,i,k,
             v,n,l,i,k,q,d,d,g,g,
             s,p,i,r,h,y,l,v,k,y,
             r,a,l,a,s,e,w,k,p,e,
             i,r,l,p,s,g,s,d,h,v,
             m,l,k,s,l,d,w,n,a,e,
             y,e,v,y,v,v,a,e,n,q,
             q,g,k,s,k,a,a,h,f,v,
             f,r,t,s,a],
  Secondary = [ strand(7,13), strand(18,24), strand(59,63), strand(51,52), 
                strand(33,42), strand(70,79), strand(82,92)].




protein(ID, Primary, Secondary):-
%% protein length: 96
  ID = '1HS7',
%% commented the short helix
  Primary = [t,n,q,k,t,k,e,l,s,n,
             l,i,e,t,f,a,e,q,s,r,
             v,l,e,k,e,c,t,k,i,g,
             s,k,r,d,s,k,e,l,r,y,
             k,i,e,t,e,l,i,p,n,c,
             t,s,v,r,d,k,i,e,s,n,
             i,l,i,h,q,n,g,k,l,s,
             a,d,f,k,n,l,k,t,k,y,
             q,s,l,q,q,s,y,n,q,r,
             k,s,l,f,p,l],
  Secondary = [ helix(2,29), helix(35,45), helix(45,60), %helix(61,66), 
               helix(66,94)].



protein(ID, Primary, Secondary):-
%% protein length: 97
  ID = '1J27',
  Primary = [m,k,a,y,l,g,l,y,t,a,
             r,l,e,t,p,a,r,s,l,k,
             e,k,r,a,l,i,k,p,a,l,
             e,r,l,k,a,r,f,p,v,s,
             a,a,r,l,y,g,l,d,a,w,
             g,y,e,v,v,g,f,t,l,l,
             g,n,d,p,a,w,v,e,e,t,
             m,r,a,a,a,r,f,l,a,e,
             a,g,g,f,q,v,a,l,e,e,
             f,r,l,e,a,f,e],
  Secondary = [ helix(18,37), helix(63,81), strand(40,45), strand(52,61), strand(2,13), strand(85,97)].



protein(ID, Primary, Secondary):-
%% protein length: 98
  ID = '2NCM',
  Primary = [r,v,l,q,v,d,i,v,p,s,
             q,g,e,i,s,v,g,e,s,k,
             f,f,l,c,q,v,a,g,d,a,
             k,d,k,d,i,s,w,f,s,p,
             n,g,e,k,l,s,p,n,q,q,
             r,i,s,v,v,w,n,d,d,d,
             s,s,t,l,t,i,y,n,a,n,
             i,d,d,a,g,i,y,k,c,v,
             v,t,a,e,d,g,t,q,s,e,
             a,t,v,n,v,k,i,f],
  Secondary = [ ssbond(24,79), strand(3,8), strand(19,28), strand(61,67), strand(52,57), strand(11,15), strand(87,98), strand(75,83), strand(34,39), strand(43,45)].


protein(ID, Primary, Secondary):-
%% protein length: 98
  ID = '1BM8',
  Primary = [q,i,y,s,a,r,y,s,g,v,
             d,v,y,e,f,i,h,s,t,g,
             s,i,m,k,r,k,k,d,d,w,
             v,n,a,t,h,i,l,k,a,a,
             n,f,a,k,a,k,r,t,r,i,
             l,e,k,e,v,l,k,e,t,h,
             e,k,v,q,g,g,f,g,k,y,
             q,g,t,w,v,p,l,n,i,a,
             k,q,l,a,e,k,f,s,v,y,
             d,q,l,k,p,l,f,d],
  Secondary = [ helix(33,39), helix(44,54), helix(77,86), helix(90,97), strand(2,7), strand(10,16), strand(21,25), strand(61,63), strand(72,74)].




protein(ID, Primary, Secondary):-
%% protein length: 100
  ID = '1JHG',
  Primary = [s,a,a,m,a,e,q,r,h,q,
             e,w,l,r,f,v,d,l,l,k,
             n,a,y,q,n,d,l,h,l,p,
             l,l,n,l,m,l,t,p,d,e,
             r,e,a,l,g,t,r,v,r,i,
             i,e,e,l,l,r,g,e,m,s,
             q,r,e,l,k,n,e,l,g,a,
             g,i,a,t,i,t,r,g,s,n,
             s,l,k,a,a,p,v,e,l,r,
             q,w,l,e,e,v,l,l,k,s],
  Secondary = [ helix(2,24), helix(28,35), helix(38,56), helix(61,68), 
                helix(72,84), helix(87,97)].



protein(ID, Primary, Secondary):-
%% protein length: 101
  ID = '1UKU',
  Primary = [m,i,i,v,y,t,t,f,p,d,
             w,e,s,a,e,k,v,v,k,t,
             l,l,k,e,r,l,i,a,c,a,
             n,l,r,e,h,r,a,f,y,w,
             w,e,g,k,i,e,e,d,k,e,
             v,g,a,i,l,k,t,r,e,d,
             l,w,e,e,l,k,e,r,i,k,
             e,l,h,p,y,d,v,p,a,i,
             i,r,i,d,v,d,d,v,n,e,
             d,y,l,k,w,l,i,e,e,t,
             k],
  Secondary = [ helix(10,24), helix(58,60), helix(61,73), helix(89,100), strand(29,41), strand(44,57), strand(2,8), strand(81,84)].


protein(ID, Primary, Secondary):-
%% protein length: 101
  ID = '1P1L',
  Primary = [m,h,n,f,i,y,i,t,a,p,
             s,l,e,e,a,e,r,i,a,k,
             r,l,l,e,k,k,l,a,a,c,
             v,n,i,f,p,i,k,s,f,f,
             w,w,e,g,k,i,e,a,a,t,
             e,f,a,m,i,v,k,t,r,s,
             e,k,f,a,e,v,r,d,e,v,
             k,a,m,h,s,y,t,t,p,c,
             i,c,a,i,p,i,e,r,g,l,
             k,e,f,l,d,w,i,d,e,t,
             v],
  Secondary = [ helix(11,25), helix(62,74), helix(90,102), strand(30,42), strand(45,59), strand(2,9), strand(81,85)].



protein(ID, Primary, Secondary):-
%% protein length: 102
  ID = '',
  Primary = [m,e,e,v,v,l,i,t,v,p,
             s,e,e,v,a,r,t,i,a,k,
             a,l,v,e,e,r,l,a,a,c,
             v,n,i,v,p,g,l,t,s,i,
             y,r,w,q,g,e,v,v,e,d,
             q,e,l,l,l,l,v,k,t,t,
             t,h,a,f,p,k,l,k,e,r,
             v,k,a,l,h,p,y,t,v,p,
             e,i,v,a,l,p,i,a,e,g,
             n,r,e,y,l,d,w,l,r,e,
             n,t],
  Secondary = [ helix(11,25), helix(63,75), helix(91,101), strand(30,42), strand(47,60), strand(2,9), strand(82,86)].


protein(ID, Primary, Secondary):-
%% protein length: 104
  ID = '1TQG',
  Primary = [g,s,h,m,e,y,l,g,v,f,
             v,d,e,t,k,e,y,l,q,n,
             l,n,d,t,l,l,e,l,e,k,
             n,p,e,d,m,e,l,i,n,e,
             a,f,r,a,l,h,t,l,k,g,
             m,a,g,t,m,g,f,s,s,m,
             a,k,l,c,h,t,l,e,n,i,
             l,d,k,a,r,n,s,e,i,k,
             i,t,s,d,l,l,d,k,i,f,
             a,g,v,d,m,i,t,r,m,v,
             d,k,i,v],
  Secondary = [ helix(2,30), helix(33,54), helix(56,75), helix(81,103)].



protein(ID, Primary, Secondary):-
%% protein length: 105
  ID = '1SA8',
  Primary = [a,f,d,g,t,w,k,v,g,g,
             l,k,l,t,i,t,q,e,g,n,
             k,f,t,v,k,e,s,s,n,f,
             r,n,i,d,v,v,f,e,l,g,
             v,d,f,a,y,s,l,a,d,g,
             t,e,l,t,g,t,w,t,m,e,
             g,n,k,l,v,g,k,f,k,r,
             v,d,n,g,k,e,l,i,a,v,
             r,e,i,s,g,n,e,l,i,q,
             t,y,t,y,e,g,v,e,a,k,
             r,i,f,k,k],
  Secondary = [ strand(4,5), strand(12,17), strand(22,28), strand(31,37), strand(41,47), strand(51,60), strand(63,70), strand(81,82), strand(77,79), strand(87,94), strand(97,104)].

protein(ID, Primary, Secondary):-
%% protein length: 109
  ID = '1TKN',
  Primary = [r,v,e,a,m,l,n,d,r,r,
             r,l,a,l,e,n,y,i,t,a,
             l,q,a,v,p,p,r,p,r,h,
             v,f,n,m,l,k,k,y,v,r,
             a,e,q,k,d,r,q,h,t,l,
             k,h,f,e,h,v,r,m,v,d,
             p,k,k,a,a,q,i,r,s,q,
             v,m,t,h,l,r,v,i,y,e,
             r,m,n,q,s,l,s,l,l,y,
             n,v,p,a,v,a,e,e,i,q,
             d,e,v,d,e,l,l,q,k],
  Secondary = [ helix(2,24), helix(27,60), helix(60,91), helix(92,108)].



protein(ID, Primary, Secondary):-
%% protein length: 116
  ID = '1NXS',
  Primary = [k,k,i,l,i,v,d,d,e,k,
             p,i,s,d,i,i,k,f,n,m,
             t,k,e,g,y,e,v,v,t,a,
             f,n,g,r,e,a,l,e,q,f,
             e,a,e,q,p,d,i,i,i,l,
             d,l,p,e,i,d,g,l,e,v,
             a,k,t,i,r,k,t,s,s,v,
             p,i,l,m,l,s,a,k,d,s,
             e,f,d,k,v,i,g,l,e,l,
             g,a,d,d,y,v,t,k,p,f,
             s,n,r,e,l,q,a,r,v,k,
             a,l,l,r,r,s],
  Secondary = [ helix(9,23), helix(32,44), helix(57,70), helix(82,92), helix(103,119), strand(26,30), strand(2,7), strand(47,51), strand(74,78), strand(96,99)].



protein(ID, Primary, Secondary):-
%% protein length: 133
  ID = '1QVX',
  Primary = [r,s,n,d,k,v,y,e,n,v,
             t,g,l,v,k,a,v,i,e,m,
             s,s,k,i,q,p,a,p,p,e,
             e,y,v,p,m,v,k,e,v,g,
             l,a,l,r,t,l,l,a,t,v,
             d,e,s,l,p,v,l,p,a,s,
             t,h,r,e,i,e,m,a,q,k,
             l,l,n,s,d,l,a,e,l,i,
             n,k,m,k,l,a,q,q,y,v,
             m,t,s,l,q,q,e,y,k,k,
             q,m,l,t,a,a,h,a,l,a,
             v,d,a,k,n,l,l,d,v,i,
             d,q,a,r,l,k,m,i,s,q,
             s,r,p],
  Secondary = [ helix(4,25), helix(31,60), helix(63,-910), helix(-906,-868)].



protein(ID, Primary, Secondary):-
%% protein length: 147
  ID = '1QJ8',
  Primary = [a,t,s,t,v,t,g,g,y,a,
             q,s,d,a,q,g,q,m,n,k,
             m,g,g,f,n,l,k,y,r,y,
             e,e,d,n,s,p,l,g,v,i,
             g,s,f,t,y,t,e,k,s,r,
             t,a,s,s,g,d,y,n,k,n,
             q,y,y,g,i,t,a,g,p,a,
             y,r,i,n,d,w,a,s,i,y,
             g,v,v,g,v,g,y,g,k,f,
             q,t,t,e,y,p,t,y,k,n,
             d,t,s,d,y,g,f,s,y,g,
             a,g,l,q,f,n,p,m,e,n,
             v,a,l,d,f,s,y,e,q,s,
             r,i,r,s,v,d,v,g,t,w,
             i,a,g,v,g,y,r],
  Secondary = [ strand(21,31), strand(2,14), strand(135,147), strand(121,132), strand(98,115), strand(77,94), strand(57,71), strand(37,51), strand(22,30)].



protein(ID, Primary, Secondary):-
%% protein length: 172
  ID = '1FEW',
  Primary = [s,l,s,s,e,a,l,m,r,r,
             a,v,s,l,v,t,d,s,t,s,
             t,f,l,s,q,t,t,y,a,l,
             i,e,a,i,t,e,y,t,k,a,
             v,y,t,l,t,s,l,y,r,q,
             y,t,s,l,l,g,k,m,n,s,
             e,e,e,d,e,v,w,q,v,i,
             i,g,a,r,a,e,m,t,s,k,
             h,q,e,y,l,k,l,e,t,t,
             w,m,t,a,v,g,l,s,e,m,
             a,a,e,a,a,y,q,t,g,a,
             d,q,a,s,i,t,a,r,n,h,
             i,q,l,v,k,l,q,v,e,e,
             v,h,q,l,s,r,k,a,e,t,
             k,l,a,e,a,q,i,e,e,l,
             k,q,k,t,q,e,e,g,e,e,
             r,a,e,s,e,q,e,a,y,l,
             r,e],
  Secondary = [ helix(2,55), helix(59,108), helix(110,171)].



protein(ID, Primary, Secondary):-
%% protein length: 180
  ID = '1DOV',
  Primary = [e,s,q,f,l,k,e,e,l,v,
             v,a,v,e,d,v,r,k,q,g,
             d,l,m,k,s,a,a,g,e,f,
             a,d,d,p,c,s,s,v,k,r,
             g,n,m,v,r,a,a,r,a,l,
             l,s,a,v,t,r,l,l,i,l,
             a,d,m,a,d,v,y,k,l,l,
             v,q,l,k,v,v,e,d,g,i,
             l,k,l,r,n,a,g,n,e,q,
             d,l,g,i,q,y,k,a,l,k,
             p,e,v,d,k,l,n,i,m,a,
             a,k,r,q,q,e,l,k,d,v,
             g,n,r,d,q,m,a,a,a,r,
             g,i,l,q,k,n,v,p,i,l,
             y,t,a,s,q,a,c,l,q,h,
             p,d,v,a,a,y,k,a,n,r,
             d,l,i,y,k,q,l,q,q,a,
             v,t,g,i,s,n,a,a,q,a],
  Secondary = [ helix(3,33), helix(36,85), helix(88,117), helix(119,150), helix(153,179)].


protein(ID, Primary, Secondary):-
%% protein length: 268
  ID = '1COV',
  Primary = [r,v,a,d,t,v,g,t,g,p,
             t,n,s,e,a,i,p,a,l,t,
             a,a,e,t,g,h,t,s,q,v,
             v,p,s,d,t,m,q,t,r,h,
             v,k,n,y,h,s,r,s,e,s,
             t,i,e,n,f,l,c,r,s,a,
             c,v,y,f,t,e,y,e,n,s,
             g,a,k,r,y,a,e,w,v,i,
             t,p,r,q,a,a,q,l,r,r,
             k,l,e,f,f,t,y,v,r,f,
             d,l,e,l,t,f,v,i,t,s,
             t,q,q,p,s,t,t,q,n,q,
             d,a,q,i,l,t,h,q,i,m,
             y,v,p,p,g,g,p,v,p,d,
             k,v,d,s,y,v,w,q,t,s,
             t,n,p,s,v,f,w,t,e,g,
             n,a,p,p,r,m,s,v,p,f,
             l,s,i,g,n,a,y,s,n,f,
             y,d,g,w,s,e,f,s,r,n,
             g,v,y,g,i,n,t,l,n,n,
             m,g,t,l,y,a,r,h,v,n,
             a,g,s,t,g,p,i,k,s,t,
             i,r,i,y,f,k,p,k,h,v,
             k,a,w,i,p,r,p,p,r,l,
             c,q,y,e,k,a,k,n,v,n,
             f,q,p,s,g,v,t,t,t,r,
             q,s,i,t,t,m,t,n],
  Secondary = [ helix(22,24), helix(32,34), helix(52,56), helix(86,94), helix(145,148), helix(195,197), strand(75,79), strand(203,208), strand(126,132), strand(154,158), strand(97,99), strand(231,233), strand(165,168), strand(101,112), strand(217,228), strand(60,68)].


protein(ID, Primary, Secondary):-
%% protein length: 28
  ID = '1DU9',
  Primary = [v,g,c,e,e,c,p,m,h,c,
             k,g,k,n,a,k,p,t,c,d,
             d,g,v,c,n,c,n,v],
  Secondary = [ ssbond(3,19), ssbond(6,24), ssbond(10,26), helix(5,10), strand(17,20), strand(23,26)].

protein(ID, Primary, Secondary):-
%% protein length: 28
  ID = '1FME',
  Primary = [e,q,y,t,a,k,y,k,g,r,
             t,f,r,n,e,k,e,l,r,d,
             f,i,e,k,f,k,g,r],
  Secondary = [ helix(14,25), strand(5,6), strand(11,12)].

protein(ID, Primary, Secondary):-
%% protein length: 27
  ID = '1E0N',
  Primary = [p,g,w,e,i,i,h,e,n,g,
             r,p,l,y,y,n,a,e,q,k,
             t,k,l,h,y,p,p],
  Secondary = [ strand(3,8), strand(11,17), strand(20,24)].

%%%%%%%%%%%%

protein(ID, Primary, Secondary):-
  ID = 'ELICA',
  Primary = [p,g,w,e,i,i,h,e,n,g],             
  Secondary = [ helix(1,10)].
  
protein(ID, Primary, Secondary):-
  ID = 'STRAND',
  Primary = [p,g,w,e,i,i,h],             
  Secondary = [ strand(1,7)].


protein(ID, Primary, Secondary):-
  ID = 'CASP',
  Primary = [m,k,i,p,k,i,y,v,e,g,e,l,n,d,g,d,r,v,a,i,e,k,d,g,n,a,i,i,f,l,e,k,d,e,e,y,s,g,n,g,k,l,l,y,q,v,i,y,d,d,l,a,k,y,m,s,l,d,t,l,k,k,d,v,l,i,q,y,p,d,k,h,t,l,t,y,l,k,a,g,t,k,l,i,s,v,p,a,e,g,y,k,v,y,p,i,m,d,f,g,f,r,v,l,k,g,y,r,l,a,t,l,e,s,k,k,g,d,l,r,y,v,n,s,p,v,s,g,t,v,i,f,m,n,e,i,p,s,e,r,a,n,y,v,f,y,m,l,e,e],
  Secondary = [].
  

