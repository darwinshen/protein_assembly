%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  The 20 amino acids with their 3 possible names  %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 12/11/2009

name_conversion("ALA", a, "alanine").
name_conversion("CYS", c, "cysteine").
name_conversion("ASP", d, "aspartic_acid").
name_conversion("GLU", e, "glutamic_acid").
name_conversion("PHE", f, "phenylalanine").
name_conversion("GLY", g, "glycine").
name_conversion("HIS", h, "histidine").
name_conversion("ILE", i, "isoleucine").
name_conversion("LYS", k, "lysine").
name_conversion("LEU", l, "leucine").
name_conversion("MET", m, "methionine").
name_conversion("ASN", n, "asparagine").
name_conversion("PRO", p, "proline").
name_conversion("GLN", q, "glutamine").
name_conversion("ARG", r, "arginine").
name_conversion("SER", s, "serine").
name_conversion("THR", t, "threonine").
name_conversion("VAL", v, "valine").
name_conversion("TRP", w, "tryptophan").
name_conversion("TYR", y, "tyrosine").
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 20x20 CONTACT POTENTIAL TABLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tablep(a,a,-1119).
tablep(a,c,-1700).
tablep(a,d,-235).
tablep(a,e,-39).
tablep(a,f,-1750).
tablep(a,g,-290).
tablep(a,h,-646).
tablep(a,i,-1872).
tablep(a,k,196).
tablep(a,l,-1728).
tablep(a,m,-1517).
tablep(a,n,-371).
tablep(a,p,-374).
tablep(a,q,-323).
tablep(a,r,-327).
tablep(a,s,-607).
tablep(a,t,-717).
tablep(a,v,-1731).
tablep(a,w,-1565).
tablep(a,y,-1318).
tablep(c,a,-1700).
tablep(c,c,-3477).
tablep(c,d,-616).
tablep(c,e,-179).
tablep(c,f,-2424).
tablep(c,g,-1101).
tablep(c,h,-1499).
tablep(c,i,-2410).
tablep(c,k,-112).
tablep(c,l,-2343).
tablep(c,m,-2240).
tablep(c,n,-788).
tablep(c,p,-1196).
tablep(c,q,-835).
tablep(c,r,-771).
tablep(c,s,-1306).
tablep(c,t,-1243).
tablep(c,v,-2258).
tablep(c,w,-2080).
tablep(c,y,-1892).
tablep(d,a,-235).
tablep(d,c,-616).
tablep(d,d,179).
tablep(d,e,634).
tablep(d,f,-482).
tablep(d,g,-97).
tablep(d,h,-664).
tablep(d,i,-402).
tablep(d,k,-176).
tablep(d,l,-291).
tablep(d,m,-409).
tablep(d,n,-344).
tablep(d,p,189).
tablep(d,q,22).
tablep(d,r,-584).
tablep(d,s,-521).
tablep(d,t,-382).
tablep(d,v,-298).
tablep(d,w,-613).
tablep(d,y,-631).
tablep(e,a,-39).
tablep(e,c,-179).
tablep(e,d,634).
tablep(e,e,933).
tablep(e,f,-419).
tablep(e,g,443).
tablep(e,h,-324).
tablep(e,i,-439).
tablep(e,k,-57).
tablep(e,l,-366).
tablep(e,m,-209).
tablep(e,n,160).
tablep(e,p,257).
tablep(e,q,179).
tablep(e,r,-374).
tablep(e,s,-161).
tablep(e,t,-192).
tablep(e,v,-335).
tablep(e,w,-624).
tablep(e,y,-453).
tablep(f,a,-1750).
tablep(f,c,-2424).
tablep(f,d,-482).
tablep(f,e,-419).
tablep(f,f,-2467).
tablep(f,g,-1034).
tablep(f,h,-1330).
tablep(f,i,-2530).
tablep(f,k,-270).
tablep(f,l,-2491).
tablep(f,m,-2304).
tablep(f,n,-790).
tablep(f,p,-1076).
tablep(f,q,-807).
tablep(f,r,-805).
tablep(f,s,-1178).
tablep(f,t,-1237).
tablep(f,v,-2391).
tablep(f,w,-2286).
tablep(f,y,-1963).
tablep(g,a,-290).
tablep(g,c,-1101).
tablep(g,d,-97).
tablep(g,e,443).
tablep(g,f,-1034).
tablep(g,g,219).
tablep(g,h,-325).
tablep(g,i,-885).
tablep(g,k,589).
tablep(g,l,-767).
tablep(g,m,-897).
tablep(g,n,-230).
tablep(g,p,-42).
tablep(g,q,33).
tablep(g,r,-50).
tablep(g,s,-261).
tablep(g,t,-311).
tablep(g,v,-756).
tablep(g,w,-1142).
tablep(g,y,-818).
tablep(h,a,-646).
tablep(h,c,-1499).
tablep(h,d,-664).
tablep(h,e,-324).
tablep(h,f,-1330).
tablep(h,g,-325).
tablep(h,h,-1078).
tablep(h,i,-1234).
tablep(h,k,388).
tablep(h,l,-1176).
tablep(h,m,-1252).
tablep(h,n,-455).
tablep(h,p,-346).
tablep(h,q,-290).
tablep(h,r,-307).
tablep(h,s,-639).
tablep(h,t,-720).
tablep(h,v,-1118).
tablep(h,w,-1383).
tablep(h,y,-1222).
tablep(i,a,-1872).
tablep(i,c,-2410).
tablep(i,d,-402).
tablep(i,e,-439).
tablep(i,f,-2530).
tablep(i,g,-885).
tablep(i,h,-1234).
tablep(i,i,-2691).
tablep(i,k,-253).
tablep(i,l,-2647).
tablep(i,m,-2286).
tablep(i,n,-669).
tablep(i,p,-991).
tablep(i,q,-778).
tablep(i,r,-854).
tablep(i,s,-1037).
tablep(i,t,-1360).
tablep(i,v,-2568).
tablep(i,w,-2303).
tablep(i,y,-1998).
tablep(k,a,196).
tablep(k,c,-112).
tablep(k,d,-176).
tablep(k,e,-57).
tablep(k,f,-270).
tablep(k,g,589).
tablep(k,h,388).
tablep(k,i,-253).
tablep(k,k,1339).
tablep(k,l,-222).
tablep(k,m,-146).
tablep(k,n,271).
tablep(k,p,661).
tablep(k,q,334).
tablep(k,r,815).
tablep(k,s,223).
tablep(k,t,155).
tablep(k,v,-200).
tablep(k,w,-391).
tablep(k,y,-349).
tablep(l,a,-1728).
tablep(l,c,-2343).
tablep(l,d,-291).
tablep(l,e,-366).
tablep(l,f,-2491).
tablep(l,g,-767).
tablep(l,h,-1176).
tablep(l,i,-2647).
tablep(l,k,-222).
tablep(l,l,-2501).
tablep(l,m,-2208).
tablep(l,n,-524).
tablep(l,p,-771).
tablep(l,q,-729).
tablep(l,r,-758).
tablep(l,s,-959).
tablep(l,t,-1202).
tablep(l,v,-2447).
tablep(l,w,-2222).
tablep(l,y,-1919).
tablep(m,a,-1517).
tablep(m,c,-2240).
tablep(m,d,-409).
tablep(m,e,-209).
tablep(m,f,-2304).
tablep(m,g,-897).
tablep(m,h,-1252).
tablep(m,i,-2286).
tablep(m,k,-146).
tablep(m,l,-2208).
tablep(m,m,-1901).
tablep(m,n,-658).
tablep(m,p,-788).
tablep(m,q,-720).
tablep(m,r,-611).
tablep(m,s,-893).
tablep(m,t,-999).
tablep(m,v,-2079).
tablep(m,w,-2090).
tablep(m,y,-1834).
tablep(n,a,-371).
tablep(n,c,-788).
tablep(n,d,-344).
tablep(n,e,160).
tablep(n,f,-790).
tablep(n,g,-230).
tablep(n,h,-455).
tablep(n,i,-669).
tablep(n,k,271).
tablep(n,l,-524).
tablep(n,m,-658).
tablep(n,n,-367).
tablep(n,p,-18).
tablep(n,q,-253).
tablep(n,r,-114).
tablep(n,s,-423).
tablep(n,t,-463).
tablep(n,v,-673).
tablep(n,w,-884).
tablep(n,y,-670).
tablep(p,a,-374).
tablep(p,c,-1196).
tablep(p,d,189).
tablep(p,e,257).
tablep(p,f,-1076).
tablep(p,g,-42).
tablep(p,h,-346).
tablep(p,i,-991).
tablep(p,k,661).
tablep(p,l,-771).
tablep(p,m,-788).
tablep(p,n,-18).
tablep(p,p,129).
tablep(p,q,-35).
tablep(p,r,-23).
tablep(p,s,-199).
tablep(p,t,-222).
tablep(p,v,-886).
tablep(p,w,-1278).
tablep(p,y,-1067).
tablep(q,a,-323).
tablep(q,c,-835).
tablep(q,d,22).
tablep(q,e,179).
tablep(q,f,-807).
tablep(q,g,33).
tablep(q,h,-290).
tablep(q,i,-778).
tablep(q,k,334).
tablep(q,l,-729).
tablep(q,m,-720).
tablep(q,n,-253).
tablep(q,p,-35).
tablep(q,q,54).
tablep(q,r,-42).
tablep(q,s,-260).
tablep(q,t,-342).
tablep(q,v,-642).
tablep(q,w,-997).
tablep(q,y,-687).
tablep(r,a,-327).
tablep(r,c,-771).
tablep(r,d,-584).
tablep(r,e,-374).
tablep(r,f,-805).
tablep(r,g,-50).
tablep(r,h,-307).
tablep(r,i,-854).
tablep(r,k,815).
tablep(r,l,-758).
tablep(r,m,-611).
tablep(r,n,-114).
tablep(r,p,-23).
tablep(r,q,-42).
tablep(r,r,200).
tablep(r,s,-264).
tablep(r,t,-247).
tablep(r,v,-664).
tablep(r,w,-912).
tablep(r,y,-745).
tablep(s,a,-607).
tablep(s,c,-1306).
tablep(s,d,-521).
tablep(s,e,-161).
tablep(s,f,-1178).
tablep(s,g,-261).
tablep(s,h,-639).
tablep(s,i,-1037).
tablep(s,k,223).
tablep(s,l,-959).
tablep(s,m,-893).
tablep(s,n,-423).
tablep(s,p,-199).
tablep(s,q,-260).
tablep(s,r,-264).
tablep(s,s,-519).
tablep(s,t,-548).
tablep(s,v,-933).
tablep(s,w,-1145).
tablep(s,y,-859).
tablep(t,a,-717).
tablep(t,c,-1243).
tablep(t,d,-382).
tablep(t,e,-192).
tablep(t,f,-1237).
tablep(t,g,-311).
tablep(t,h,-720).
tablep(t,i,-1360).
tablep(t,k,155).
tablep(t,l,-1202).
tablep(t,m,-999).
tablep(t,n,-463).
tablep(t,p,-222).
tablep(t,q,-342).
tablep(t,r,-247).
tablep(t,s,-548).
tablep(t,t,-617).
tablep(t,v,-1240).
tablep(t,w,-1077).
tablep(t,y,-892).
tablep(v,a,-1731).
tablep(v,c,-2258).
tablep(v,d,-298).
tablep(v,e,-335).
tablep(v,f,-2391).
tablep(v,g,-756).
tablep(v,h,-1118).
tablep(v,i,-2568).
tablep(v,k,-200).
tablep(v,l,-2447).
tablep(v,m,-2079).
tablep(v,n,-673).
tablep(v,p,-886).
tablep(v,q,-642).
tablep(v,r,-664).
tablep(v,s,-933).
tablep(v,t,-1240).
tablep(v,v,-2385).
tablep(v,w,-2097).
tablep(v,y,-1790).
tablep(w,a,-1565).
tablep(w,c,-2080).
tablep(w,d,-613).
tablep(w,e,-624).
tablep(w,f,-2286).
tablep(w,g,-1142).
tablep(w,h,-1383).
tablep(w,i,-2303).
tablep(w,k,-391).
tablep(w,l,-2222).
tablep(w,m,-2090).
tablep(w,n,-884).
tablep(w,p,-1278).
tablep(w,q,-997).
tablep(w,r,-912).
tablep(w,s,-1145).
tablep(w,t,-1077).
tablep(w,v,-2097).
tablep(w,w,-1867).
tablep(w,y,-1834).
tablep(y,a,-1318).
tablep(y,c,-1892).
tablep(y,d,-631).
tablep(y,e,-453).
tablep(y,f,-1963).
tablep(y,g,-818).
tablep(y,h,-1222).
tablep(y,i,-1998).
tablep(y,k,-349).
tablep(y,l,-1919).
tablep(y,m,-1834).
tablep(y,n,-670).
tablep(y,p,-1067).
tablep(y,q,-687).
tablep(y,r,-745).
tablep(y,s,-859).
tablep(y,t,-892).
tablep(y,v,-1790).
tablep(y,w,-1834).
tablep(y,y,-1335).

%%%%%%%%%%%%%%%%%%%%%%
%%% Centroid radius
%%%%%%%%%%%%%%%%%%%%%%

radius(a  ,  2.40).
radius(c  ,  3.23).
radius(d  ,  3.40).
radius(e  ,  4.03).
radius(f  ,  4.26).
radius(g  ,  2.40).
radius(h  ,  4.04).
radius(i  ,  3.16).
radius(k  ,  4.41).
radius(l  ,  3.48).
radius(m  ,  4.12).
radius(n  ,  3.40).
radius(p  ,  2.70).
radius(q  ,  3.98).
radius(s  ,  2.80).
radius(t  ,  2.80).
radius(r  ,  5.01).
radius(v  ,  2.80).
radius(w  ,  4.83).
radius(y  ,  4.73).

%%%%%%%%%%%%%%%%%%%%%%
%%% Centroid Angles
%%%%%%%%%%%%%%%%%%%%%%

centroid_ang1(a, 120.92).
centroid_ang1(r, 114.83).
centroid_ang1(n, 115.73).
centroid_ang1(d, 117.74).
centroid_ang1(c, 113.60).
centroid_ang1(q, 115.21).
centroid_ang1(e, 115.80).
centroid_ang1(g, 0).
centroid_ang1(h, 112.88).
centroid_ang1(i, 119.48).
centroid_ang1(l, 120.39).
centroid_ang1(k, 116.74).
centroid_ang1(m, 116.42).
centroid_ang1(f, 112.99).
centroid_ang1(p, 81.30).
centroid_ang1(s, 116.84).
centroid_ang1(t, 115.53).
centroid_ang1(w, 113.95).
centroid_ang1(y, 112.01).
centroid_ang1(v, 119.62).

centroid_ang2(a, 110.53).
centroid_ang2(r, 113.59).
centroid_ang2(n, 117.73).
centroid_ang2(d, 116.03).
centroid_ang2(c, 115.36).
centroid_ang2(q, 115.96).
centroid_ang2(e, 115.98).
centroid_ang2(g, 0).
centroid_ang2(h, 115.38).
centroid_ang2(i, 118.17).
centroid_ang2(l, 119.90).
centroid_ang2(k, 115.73).
centroid_ang2(m, 115.79).
centroid_ang2(f, 114.40).
centroid_ang2(p, 123.58).
centroid_ang2(s, 110.33).
centroid_ang2(t, 111.67).
centroid_ang2(w, 109.27).
centroid_ang2(y, 113.14).
centroid_ang2(v, 114.46).

centroid_tors(a, -138.45).
centroid_tors(r, -155.07).
centroid_tors(n, -144.56).
centroid_tors(d, -142.28).
centroid_tors(c, -142.28).
centroid_tors(q, -149.99).
centroid_tors(e, -147.56).
centroid_tors(g, -0     ).
centroid_tors(h, -144.08).
centroid_tors(i, -151.72).
centroid_tors(l, -153.24).
centroid_tors(k, -153.03).
centroid_tors(m, -159.50).
centroid_tors(f, -146.92).
centroid_tors(p, -105.62).
centroid_tors(s, -139.94).
centroid_tors(t, -142.28).
centroid_tors(w, -155.35).
centroid_tors(y, -149.29).
centroid_tors(v, -150.47).

centroid_dist(a, 1.55).
centroid_dist(r, 4.16).
centroid_dist(n, 2.55).
centroid_dist(d, 2.55).
centroid_dist(c, 2.38).
centroid_dist(q, 3.13).
centroid_dist(e, 3.18).
centroid_dist(g, 0  ).
centroid_dist(h, 3.19).
centroid_dist(i, 2.31).
centroid_dist(l, 2.63).
centroid_dist(k, 3.56).
centroid_dist(m, 3.27).
centroid_dist(f, 3.41).
centroid_dist(p, 1.85).
centroid_dist(s, 1.95).
centroid_dist(t, 1.95).
centroid_dist(w, 3.98).
centroid_dist(y, 3.88).
centroid_dist(v, 1.95).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%     END        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

