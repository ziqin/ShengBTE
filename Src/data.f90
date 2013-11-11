!  ShengBTE, a solver for the Boltzmann Transport Equation for phonons
!  Copyright (C) 2012-2013 Wu Li <wu.li@cea.fr>
!  Copyright (C) 2012-2013 Jesús Carrete Montaña <jcarrete@gmail.com>
!  Copyright (C) 2012-2013 Nebil Ayape Katcho <nebil.ayapekatcho@cea.fr>
!  Copyright (C) 2012-2013 Natalio Mingo Bisquert <natalio.mingo@cea.fr>
!
!  This program is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.

! Physcal constants and other relevant data.
module data
  implicit none
  
  real(kind=8),parameter :: pi=3.141592653589793238d0
  real(kind=8),parameter :: kb=1.380648813d-23 ! J/K
  real(kind=8),parameter :: hbar=1.05457172647d-22 ! J*THz
  complex(kind=8),parameter :: iunit=(0.,1.)

  character(len=3),parameter :: periodic_table(114)=[character(len=3) ::&
       "H","He","Li","Be","B","C","N","O",&
       "F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V",&
       "Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb",&
       "Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb",&
       "Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb",&
       "Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",&
       "Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am",&
       "Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds",&
       "Rg","Cn","Uuq","Uuh"]
  
  type :: isotope_t
     character(len=3) :: element
     real(kind=8) :: mass
     real(kind=8) :: abundance
     type(isotope_t),pointer :: prev=>null()
  end type isotope_t

  type(isotope_t),pointer :: isotope_end=>null()

contains

  subroutine data_new_isotope(element,mass,abundance)
    character(len=*),intent(in) :: element
    real(kind=8),intent(in) :: mass,abundance
    type(isotope_t),pointer :: new_isotope


    allocate(new_isotope)
    new_isotope%mass=mass
    new_isotope%abundance=abundance
    new_isotope%element=element
    new_isotope%prev=>isotope_end
    isotope_end=>new_isotope
  end subroutine data_new_isotope
    
  subroutine data_fill_isotopes()
    call data_new_isotope("Ag",106.905095d0,51.84d0)
    call data_new_isotope("Ag",108.904754d0,48.16d0)
    call data_new_isotope("Al",26.981541d0,100.0d0)
    call data_new_isotope("Ar",35.967546d0,0.34d0)
    call data_new_isotope("Ar",37.962732d0,0.063d0)
    call data_new_isotope("Ar",39.962383d0,99.6d0)
    call data_new_isotope("As",74.921596d0,100.0d0)
    call data_new_isotope("Au",196.96656d0,100.0d0)
    call data_new_isotope("B",10.012938d0,19.8d0)
    call data_new_isotope("B",11.009305d0,80.2d0)
    call data_new_isotope("Ba",129.906277d0,0.11d0)
    call data_new_isotope("Ba",131.905042d0,0.1d0)
    call data_new_isotope("Ba",133.90449d0,2.42d0)
    call data_new_isotope("Ba",134.905668d0,6.59d0)
    call data_new_isotope("Ba",135.904556d0,7.85d0)
    call data_new_isotope("Ba",136.905816d0,11.23d0)
    call data_new_isotope("Ba",137.905236d0,71.7d0)
    call data_new_isotope("Be",9.012183d0,100.0d0)
    call data_new_isotope("Bi",208.980388d0,100.0d0)
    call data_new_isotope("Br",78.918336d0,50.69d0)
    call data_new_isotope("Br",80.91629d0,49.31d0)
    call data_new_isotope("C",12.0d0,98.9d0)
    call data_new_isotope("C",13.003355d0,1.1d0)
    call data_new_isotope("Ca",39.962591d0,96.95d0)
    call data_new_isotope("Ca",41.958622d0,0.65d0)
    call data_new_isotope("Ca",42.95877d0,0.14d0)
    call data_new_isotope("Ca",43.955485d0,2.086d0)
    call data_new_isotope("Ca",45.953689d0,0.004d0)
    call data_new_isotope("Ca",47.952532d0,0.19d0)
    call data_new_isotope("Cd",105.906461d0,1.25d0)
    call data_new_isotope("Cd",107.904186d0,0.89d0)
    call data_new_isotope("Cd",109.903007d0,12.49d0)
    call data_new_isotope("Cd",110.904182d0,12.8d0)
    call data_new_isotope("Cd",111.902761d0,24.13d0)
    call data_new_isotope("Cd",112.904401d0,12.22d0)
    call data_new_isotope("Cd",113.903361d0,28.73d0)
    call data_new_isotope("Cd",115.904758d0,7.49d0)
    call data_new_isotope("Ce",135.90714d0,0.19d0)
    call data_new_isotope("Ce",137.905996d0,0.25d0)
    call data_new_isotope("Ce",139.905442d0,88.48d0)
    call data_new_isotope("Ce",141.909249d0,11.08d0)
    call data_new_isotope("Cl",34.968853d0,75.77d0)
    call data_new_isotope("Cl",36.965903d0,24.23d0)
    call data_new_isotope("Co",58.933198d0,100.0d0)
    call data_new_isotope("Cr",49.946046d0,4.35d0)
    call data_new_isotope("Cr",51.94051d0,83.79d0)
    call data_new_isotope("Cr",52.940651d0,9.5d0)
    call data_new_isotope("Cr",53.938882d0,2.36d0)
    call data_new_isotope("Cs",132.905433d0,100.0d0)
    call data_new_isotope("Cu",62.929599d0,69.17d0)
    call data_new_isotope("Cu",64.927792d0,30.83d0)
    call data_new_isotope("Dy",155.924287d0,0.06d0)
    call data_new_isotope("Dy",157.924412d0,0.1d0)
    call data_new_isotope("Dy",159.925203d0,2.34d0)
    call data_new_isotope("Dy",160.926939d0,18.9d0)
    call data_new_isotope("Dy",161.926805d0,25.5d0)
    call data_new_isotope("Dy",162.928737d0,24.9d0)
    call data_new_isotope("Dy",163.929183d0,28.2d0)
    call data_new_isotope("Er",161.928787d0,0.14d0)
    call data_new_isotope("Er",163.929211d0,1.61d0)
    call data_new_isotope("Er",165.930305d0,33.6d0)
    call data_new_isotope("Er",166.932061d0,22.95d0)
    call data_new_isotope("Er",167.932383d0,26.8d0)
    call data_new_isotope("Er",169.935476d0,14.9d0)
    call data_new_isotope("Eu",150.91986d0,47.8d0)
    call data_new_isotope("Eu",152.921243d0,52.2d0)
    call data_new_isotope("F",18.998403d0,100.0d0)
    call data_new_isotope("Fe",53.939612d0,5.8d0)
    call data_new_isotope("Fe",55.934939d0,91.72d0)
    call data_new_isotope("Fe",56.935396d0,2.2d0)
    call data_new_isotope("Fe",57.933278d0,0.28d0)
    call data_new_isotope("Ga",68.925581d0,60.1d0)
    call data_new_isotope("Ga",70.924701d0,39.9d0)
    call data_new_isotope("Gd",151.919803d0,0.2d0)
    call data_new_isotope("Gd",153.920876d0,2.18d0)
    call data_new_isotope("Gd",154.822629d0,14.8d0)
    call data_new_isotope("Gd",155.92213d0,20.47d0)
    call data_new_isotope("Gd",156.923967d0,15.65d0)
    call data_new_isotope("Gd",157.924111d0,24.84d0)
    call data_new_isotope("Gd",159.927061d0,21.86d0)
    call data_new_isotope("Ge",69.92425d0,20.5d0)
    call data_new_isotope("Ge",71.92208d0,27.4d0)
    call data_new_isotope("Ge",72.923464d0,7.8d0)
    call data_new_isotope("Ge",73.921179d0,36.5d0)
    call data_new_isotope("Ge",75.921403d0,7.8d0)
    call data_new_isotope("H",1.007825d0,99.99d0)
    call data_new_isotope("H",2.014102d0,0.015d0)
    call data_new_isotope("He",3.016029d0,0.0001d0)
    call data_new_isotope("He",4.002603d0,100.0d0)
    call data_new_isotope("Hf",173.940065d0,0.16d0)
    call data_new_isotope("Hf",175.94142d0,5.2d0)
    call data_new_isotope("Hf",176.943233d0,18.6d0)
    call data_new_isotope("Hf",177.94371d0,27.1d0)
    call data_new_isotope("Hf",178.945827d0,13.74d0)
    call data_new_isotope("Hf",179.946561d0,35.2d0)
    call data_new_isotope("Hg",195.965812d0,0.15d0)
    call data_new_isotope("Hg",197.96676d0,10.1d0)
    call data_new_isotope("Hg",198.968269d0,17.0d0)
    call data_new_isotope("Hg",199.968316d0,23.1d0)
    call data_new_isotope("Hg",200.970293d0,13.2d0)
    call data_new_isotope("Hg",201.970632d0,29.65d0)
    call data_new_isotope("Hg",203.973481d0,6.8d0)
    call data_new_isotope("Ho",164.930332d0,100.0d0)
    call data_new_isotope("I",126.904477d0,100.0d0)
    call data_new_isotope("In",112.904056d0,4.3d0)
    call data_new_isotope("In",114.903875d0,95.7d0)
    call data_new_isotope("Ir",190.960603d0,37.3d0)
    call data_new_isotope("Ir",192.962942d0,62.7d0)
    call data_new_isotope("K",38.963708d0,93.2d0)
    call data_new_isotope("K",39.963999d0,0.012d0)
    call data_new_isotope("K",40.961825d0,6.73d0)
    call data_new_isotope("Kr",77.920397d0,0.35d0)
    call data_new_isotope("Kr",79.916375d0,2.25d0)
    call data_new_isotope("Kr",81.913483d0,11.6d0)
    call data_new_isotope("Kr",82.914134d0,11.5d0)
    call data_new_isotope("Kr",83.911506d0,57.0d0)
    call data_new_isotope("Kr",85.910614d0,17.3d0)
    call data_new_isotope("La",137.907114d0,0.09d0)
    call data_new_isotope("La",138.906355d0,99.91d0)
    call data_new_isotope("Li",6.015123d0,7.42d0)
    call data_new_isotope("Li",7.016005d0,92.58d0)
    call data_new_isotope("Lu",174.940785d0,97.4d0)
    call data_new_isotope("Lu",175.942694d0,2.6d0)
    call data_new_isotope("Mg",23.985045d0,78.9d0)
    call data_new_isotope("Mg",24.985839d0,10.0d0)
    call data_new_isotope("Mg",25.982595d0,11.1d0)
    call data_new_isotope("Mn",54.938046d0,100.0d0)
    call data_new_isotope("Mo",91.906809d0,14.84d0)
    call data_new_isotope("Mo",93.905086d0,9.25d0)
    call data_new_isotope("Mo",94.905838d0,15.92d0)
    call data_new_isotope("Mo",95.904676d0,16.68d0)
    call data_new_isotope("Mo",96.906018d0,9.55d0)
    call data_new_isotope("Mo",97.905405d0,24.13d0)
    call data_new_isotope("Mo",99.907473d0,9.63d0)
    call data_new_isotope("N",14.003074d0,99.63d0)
    call data_new_isotope("N",15.000109d0,0.37d0)
    call data_new_isotope("Na",22.98977d0,100.0d0)
    call data_new_isotope("Nb",92.906378d0,100.0d0)
    call data_new_isotope("Nd",141.907731d0,27.13d0)
    call data_new_isotope("Nd",142.909823d0,12.18d0)
    call data_new_isotope("Nd",143.910096d0,23.8d0)
    call data_new_isotope("Nd",144.912582d0,8.3d0)
    call data_new_isotope("Nd",145.913126d0,17.19d0)
    call data_new_isotope("Nd",147.916901d0,5.76d0)
    call data_new_isotope("Nd",149.9209d0,5.64d0)
    call data_new_isotope("Ne",19.992439d0,90.6d0)
    call data_new_isotope("Ne",20.993845d0,0.26d0)
    call data_new_isotope("Ne",21.991384d0,9.2d0)
    call data_new_isotope("Ni",57.935347d0,68.27d0)
    call data_new_isotope("Ni",59.930789d0,26.1d0)
    call data_new_isotope("Ni",60.931059d0,1.13d0)
    call data_new_isotope("Ni",61.928346d0,3.59d0)
    call data_new_isotope("Ni",63.927968d0,0.91d0)
    call data_new_isotope("O",15.994915d0,99.76d0)
    call data_new_isotope("O",16.999131d0,0.038d0)
    call data_new_isotope("O",17.999159d0,0.2d0)
    call data_new_isotope("Os",183.952514d0,0.02d0)
    call data_new_isotope("Os",185.953852d0,1.58d0)
    call data_new_isotope("Os",186.955762d0,1.6d0)
    call data_new_isotope("Os",187.95585d0,13.3d0)
    call data_new_isotope("Os",188.958156d0,16.1d0)
    call data_new_isotope("Os",189.958455d0,26.4d0)
    call data_new_isotope("Os",191.961487d0,41.0d0)
    call data_new_isotope("P",30.973763d0,100.0d0)
    call data_new_isotope("Pb",203.973037d0,1.4d0)
    call data_new_isotope("Pb",205.974455d0,24.1d0)
    call data_new_isotope("Pb",206.975885d0,22.1d0)
    call data_new_isotope("Pb",207.976641d0,52.4d0)
    call data_new_isotope("Pd",101.905609d0,1.02d0)
    call data_new_isotope("Pd",103.904026d0,11.14d0)
    call data_new_isotope("Pd",104.905075d0,22.33d0)
    call data_new_isotope("Pd",105.903475d0,27.33d0)
    call data_new_isotope("Pd",107.903894d0,26.46d0)
    call data_new_isotope("Pd",109.905169d0,11.72d0)
    call data_new_isotope("Pr",140.907657d0,100.0d0)
    call data_new_isotope("Pt",189.959937d0,0.01d0)
    call data_new_isotope("Pt",191.961049d0,0.79d0)
    call data_new_isotope("Pt",193.962679d0,32.9d0)
    call data_new_isotope("Pt",194.964785d0,33.8d0)
    call data_new_isotope("Pt",195.964947d0,25.3d0)
    call data_new_isotope("Pt",197.967879d0,7.2d0)
    call data_new_isotope("Rb",84.9118d0,72.17d0)
    call data_new_isotope("Rb",86.909184d0,27.84d0)
    call data_new_isotope("Re",184.952977d0,37.4d0)
    call data_new_isotope("Re",186.955765d0,62.6d0)
    call data_new_isotope("Rh",102.905503d0,100.0d0)
    call data_new_isotope("Ru",95.907596d0,5.52d0)
    call data_new_isotope("Ru",97.905287d0,1.88d0)
    call data_new_isotope("Ru",98.905937d0,12.7d0)
    call data_new_isotope("Ru",99.904218d0,12.6d0)
    call data_new_isotope("Ru",100.905581d0,17.0d0)
    call data_new_isotope("Ru",101.904348d0,31.6d0)
    call data_new_isotope("Ru",103.905422d0,18.7d0)
    call data_new_isotope("S",31.972072d0,95.02d0)
    call data_new_isotope("S",32.971459d0,0.75d0)
    call data_new_isotope("S",33.967868d0,4.21d0)
    call data_new_isotope("S",35.967079d0,0.02d0)
    call data_new_isotope("Sb",120.903824d0,57.3d0)
    call data_new_isotope("Sb",122.904222d0,42.7d0)
    call data_new_isotope("Sc",44.955914d0,100.0d0)
    call data_new_isotope("Se",73.922477d0,0.9d0)
    call data_new_isotope("Se",75.919207d0,9.0d0)
    call data_new_isotope("Se",76.919908d0,7.6d0)
    call data_new_isotope("Se",77.917304d0,23.5d0)
    call data_new_isotope("Se",79.916521d0,49.6d0)
    call data_new_isotope("Se",81.916709d0,9.4d0)
    call data_new_isotope("Si",27.976928d0,92.23d0)
    call data_new_isotope("Si",28.976496d0,4.67d0)
    call data_new_isotope("Si",29.973772d0,3.1d0)
    call data_new_isotope("Sm",143.912009d0,3.1d0)
    call data_new_isotope("Sm",146.914907d0,15.0d0)
    call data_new_isotope("Sm",147.914832d0,11.3d0)
    call data_new_isotope("Sm",148.917193d0,13.8d0)
    call data_new_isotope("Sm",149.917285d0,7.4d0)
    call data_new_isotope("Sm",151.919741d0,26.7d0)
    call data_new_isotope("Sm",153.922218d0,22.7d0)
    call data_new_isotope("Sn",111.904826d0,0.97d0)
    call data_new_isotope("Sn",113.902784d0,0.65d0)
    call data_new_isotope("Sn",114.903348d0,0.36d0)
    call data_new_isotope("Sn",115.901744d0,14.7d0)
    call data_new_isotope("Sn",116.902954d0,7.7d0)
    call data_new_isotope("Sn",117.901607d0,24.3d0)
    call data_new_isotope("Sn",118.90331d0,8.6d0)
    call data_new_isotope("Sn",119.902199d0,32.4d0)
    call data_new_isotope("Sn",121.90344d0,4.6d0)
    call data_new_isotope("Sn",123.905271d0,5.6d0)
    call data_new_isotope("Sr",83.913428d0,0.56d0)
    call data_new_isotope("Sr",85.909273d0,9.86d0)
    call data_new_isotope("Sr",86.908902d0,7.0d0)
    call data_new_isotope("Sr",87.905625d0,82.58d0)
    call data_new_isotope("Ta",179.947489d0,0.012d0)
    call data_new_isotope("Ta",180.948014d0,99.99d0)
    call data_new_isotope("Tb",158.92535d0,100.0d0)
    call data_new_isotope("Te",119.904021d0,0.096d0)
    call data_new_isotope("Te",121.903055d0,2.6d0)
    call data_new_isotope("Te",122.904278d0,0.91d0)
    call data_new_isotope("Te",123.902825d0,4.82d0)
    call data_new_isotope("Te",124.904435d0,7.14d0)
    call data_new_isotope("Te",125.90331d0,18.95d0)
    call data_new_isotope("Te",127.904464d0,31.69d0)
    call data_new_isotope("Te",129.906229d0,33.8d0)
    call data_new_isotope("Th",232.038054d0,100.0d0)
    call data_new_isotope("Ti",45.952633d0,8.0d0)
    call data_new_isotope("Ti",46.951765d0,7.3d0)
    call data_new_isotope("Ti",47.947947d0,73.8d0)
    call data_new_isotope("Ti",48.947871d0,5.5d0)
    call data_new_isotope("Ti",49.944786d0,5.4d0)
    call data_new_isotope("Tl",202.972336d0,29.52d0)
    call data_new_isotope("Tl",204.97441d0,70.48d0)
    call data_new_isotope("Tm",168.934225d0,100.0d0)
    call data_new_isotope("U",234.040947d0,0.006d0)
    call data_new_isotope("U",235.043925d0,0.72d0)
    call data_new_isotope("U",238.050786d0,99.27d0)
    call data_new_isotope("V",49.947161d0,0.25d0)
    call data_new_isotope("V",50.943963d0,99.75d0)
    call data_new_isotope("W",179.946727d0,0.13d0)
    call data_new_isotope("W",181.948225d0,26.3d0)
    call data_new_isotope("W",182.950245d0,14.3d0)
    call data_new_isotope("W",183.950953d0,30.67d0)
    call data_new_isotope("W",185.954377d0,28.6d0)
    call data_new_isotope("Xe",123.905894d0,0.1d0)
    call data_new_isotope("Xe",125.904281d0,0.09d0)
    call data_new_isotope("Xe",127.903531d0,1.91d0)
    call data_new_isotope("Xe",128.90478d0,26.4d0)
    call data_new_isotope("Xe",129.90351d0,4.1d0)
    call data_new_isotope("Xe",130.905076d0,21.2d0)
    call data_new_isotope("Xe",131.904148d0,26.9d0)
    call data_new_isotope("Xe",133.905395d0,10.4d0)
    call data_new_isotope("Xe",135.907219d0,8.9d0)
    call data_new_isotope("Y",88.905856d0,100.0d0)
    call data_new_isotope("Yb",167.933908d0,0.13d0)
    call data_new_isotope("Yb",169.934774d0,3.05d0)
    call data_new_isotope("Yb",170.936338d0,14.3d0)
    call data_new_isotope("Yb",171.936393d0,21.9d0)
    call data_new_isotope("Yb",172.938222d0,16.12d0)
    call data_new_isotope("Yb",173.938873d0,31.8d0)
    call data_new_isotope("Yb",175.942576d0,12.7d0)
    call data_new_isotope("Zn",63.929145d0,48.6d0)
    call data_new_isotope("Zn",65.926035d0,27.9d0)
    call data_new_isotope("Zn",66.927129d0,4.1d0)
    call data_new_isotope("Zn",67.924846d0,18.8d0)
    call data_new_isotope("Zn",69.925325d0,0.6d0)
    call data_new_isotope("Zr",89.904708d0,51.45d0)
    call data_new_isotope("Zr",90.905644d0,11.27d0)
    call data_new_isotope("Zr",91.905039d0,17.17d0)
    call data_new_isotope("Zr",93.906319d0,17.33d0)
    call data_new_isotope("Zr",95.908272d0,2.78d0)
  end subroutine data_fill_isotopes

  subroutine data_free_isotopes()
    type(isotope_t),pointer :: p

    do
       if(.not.associated(isotope_end))exit
       p=>isotope_end%prev
       deallocate(isotope_end)
       isotope_end=>p
    end do
  end subroutine data_free_isotopes

  subroutine data_calc_mandg(element,m,g)
    character(len=3),intent(in) :: element
    real(kind=8),intent(out) :: m,g

    type(isotope_t),pointer :: p

    m=0.
    p=>isotope_end
    do
       if(.not.associated(p))exit
       if(p%element.eq.element) then
          m=m+p%mass*p%abundance
       end if
       p=>p%prev
    end do
    m=m/100.
    g=0.
    p=>isotope_end
    do
       if(.not.associated(p))exit
       if(p%element.eq.element) then
          g=g+(1.-p%mass/m)**2*p%abundance
       end if
       p=>p%prev
    end do
    g=g/100.
  end subroutine data_calc_mandg
end module data
