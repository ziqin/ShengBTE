!  ShengBTE, a solver for the Boltzmann Transport Equation for phonons
!  Copyright (C) 2012-2015 Wu Li <wu.li.phys2011@gmail.com>
!  Copyright (C) 2012-2015 Jesús Carrete Montaña <jcarrete@gmail.com>
!  Copyright (C) 2012-2015 Nebil Ayape Katcho <nebil.ayapekatcho@cea.fr>
!  Copyright (C) 2012-2015 Natalio Mingo Bisquert <natalio.mingo@cea.fr>
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

! Physical constants, table of isotopic masses and other relevant data.
module data
  implicit none

  real(kind=8),parameter :: pi=3.141592653589793238d0
  real(kind=8),parameter :: kb=1.380648813d-23 ! J/K
  real(kind=8),parameter :: hbar=1.05457172647d-22 ! J*THz
  complex(kind=8),parameter :: iunit=(0.,1.)

  ! List of elements, ordered by atomic number.
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

  ! Information about isotopes.
  character(len=3),allocatable :: isotope_element(:)
  real(kind=8),allocatable :: isotope_mass(:)
  real(kind=8),allocatable :: isotope_abundance(:)

contains

  ! Fill in the information about isotopes. This code
  ! was generated automatically.
  subroutine data_fill_isotopes()
    allocate(isotope_element(287))
    allocate(isotope_mass(287))
    allocate(isotope_abundance(287))
    isotope_element(1)="Ag"
    isotope_mass(1)=106.905095d0
    isotope_abundance(1)=51.84d0
    isotope_element(2)="Ag"
    isotope_mass(2)=108.904754d0
    isotope_abundance(2)=48.16d0
    isotope_element(3)="Al"
    isotope_mass(3)=26.981541d0
    isotope_abundance(3)=100.0d0
    isotope_element(4)="Ar"
    isotope_mass(4)=35.967546d0
    isotope_abundance(4)=0.34d0
    isotope_element(5)="Ar"
    isotope_mass(5)=37.962732d0
    isotope_abundance(5)=0.063d0
    isotope_element(6)="Ar"
    isotope_mass(6)=39.962383d0
    isotope_abundance(6)=99.6d0
    isotope_element(7)="As"
    isotope_mass(7)=74.921596d0
    isotope_abundance(7)=100.0d0
    isotope_element(8)="Au"
    isotope_mass(8)=196.96656d0
    isotope_abundance(8)=100.0d0
    isotope_element(9)="B"
    isotope_mass(9)=10.012938d0
    isotope_abundance(9)=19.8d0
    isotope_element(10)="B"
    isotope_mass(10)=11.009305d0
    isotope_abundance(10)=80.2d0
    isotope_element(11)="Ba"
    isotope_mass(11)=129.906277d0
    isotope_abundance(11)=0.11d0
    isotope_element(12)="Ba"
    isotope_mass(12)=131.905042d0
    isotope_abundance(12)=0.1d0
    isotope_element(13)="Ba"
    isotope_mass(13)=133.90449d0
    isotope_abundance(13)=2.42d0
    isotope_element(14)="Ba"
    isotope_mass(14)=134.905668d0
    isotope_abundance(14)=6.59d0
    isotope_element(15)="Ba"
    isotope_mass(15)=135.904556d0
    isotope_abundance(15)=7.85d0
    isotope_element(16)="Ba"
    isotope_mass(16)=136.905816d0
    isotope_abundance(16)=11.23d0
    isotope_element(17)="Ba"
    isotope_mass(17)=137.905236d0
    isotope_abundance(17)=71.7d0
    isotope_element(18)="Be"
    isotope_mass(18)=9.012183d0
    isotope_abundance(18)=100.0d0
    isotope_element(19)="Bi"
    isotope_mass(19)=208.980388d0
    isotope_abundance(19)=100.0d0
    isotope_element(20)="Br"
    isotope_mass(20)=78.918336d0
    isotope_abundance(20)=50.69d0
    isotope_element(21)="Br"
    isotope_mass(21)=80.91629d0
    isotope_abundance(21)=49.31d0
    isotope_element(22)="C"
    isotope_mass(22)=12.0d0
    isotope_abundance(22)=98.9d0
    isotope_element(23)="C"
    isotope_mass(23)=13.003355d0
    isotope_abundance(23)=1.1d0
    isotope_element(24)="Ca"
    isotope_mass(24)=39.962591d0
    isotope_abundance(24)=96.95d0
    isotope_element(25)="Ca"
    isotope_mass(25)=41.958622d0
    isotope_abundance(25)=0.65d0
    isotope_element(26)="Ca"
    isotope_mass(26)=42.95877d0
    isotope_abundance(26)=0.14d0
    isotope_element(27)="Ca"
    isotope_mass(27)=43.955485d0
    isotope_abundance(27)=2.086d0
    isotope_element(28)="Ca"
    isotope_mass(28)=45.953689d0
    isotope_abundance(28)=0.004d0
    isotope_element(29)="Ca"
    isotope_mass(29)=47.952532d0
    isotope_abundance(29)=0.19d0
    isotope_element(30)="Cd"
    isotope_mass(30)=105.906461d0
    isotope_abundance(30)=1.25d0
    isotope_element(31)="Cd"
    isotope_mass(31)=107.904186d0
    isotope_abundance(31)=0.89d0
    isotope_element(32)="Cd"
    isotope_mass(32)=109.903007d0
    isotope_abundance(32)=12.49d0
    isotope_element(33)="Cd"
    isotope_mass(33)=110.904182d0
    isotope_abundance(33)=12.8d0
    isotope_element(34)="Cd"
    isotope_mass(34)=111.902761d0
    isotope_abundance(34)=24.13d0
    isotope_element(35)="Cd"
    isotope_mass(35)=112.904401d0
    isotope_abundance(35)=12.22d0
    isotope_element(36)="Cd"
    isotope_mass(36)=113.903361d0
    isotope_abundance(36)=28.73d0
    isotope_element(37)="Cd"
    isotope_mass(37)=115.904758d0
    isotope_abundance(37)=7.49d0
    isotope_element(38)="Ce"
    isotope_mass(38)=135.90714d0
    isotope_abundance(38)=0.19d0
    isotope_element(39)="Ce"
    isotope_mass(39)=137.905996d0
    isotope_abundance(39)=0.25d0
    isotope_element(40)="Ce"
    isotope_mass(40)=139.905442d0
    isotope_abundance(40)=88.48d0
    isotope_element(41)="Ce"
    isotope_mass(41)=141.909249d0
    isotope_abundance(41)=11.08d0
    isotope_element(42)="Cl"
    isotope_mass(42)=34.968853d0
    isotope_abundance(42)=75.77d0
    isotope_element(43)="Cl"
    isotope_mass(43)=36.965903d0
    isotope_abundance(43)=24.23d0
    isotope_element(44)="Co"
    isotope_mass(44)=58.933198d0
    isotope_abundance(44)=100.0d0
    isotope_element(45)="Cr"
    isotope_mass(45)=49.946046d0
    isotope_abundance(45)=4.35d0
    isotope_element(46)="Cr"
    isotope_mass(46)=51.94051d0
    isotope_abundance(46)=83.79d0
    isotope_element(47)="Cr"
    isotope_mass(47)=52.940651d0
    isotope_abundance(47)=9.5d0
    isotope_element(48)="Cr"
    isotope_mass(48)=53.938882d0
    isotope_abundance(48)=2.36d0
    isotope_element(49)="Cs"
    isotope_mass(49)=132.905433d0
    isotope_abundance(49)=100.0d0
    isotope_element(50)="Cu"
    isotope_mass(50)=62.929599d0
    isotope_abundance(50)=69.17d0
    isotope_element(51)="Cu"
    isotope_mass(51)=64.927792d0
    isotope_abundance(51)=30.83d0
    isotope_element(52)="Dy"
    isotope_mass(52)=155.924287d0
    isotope_abundance(52)=0.06d0
    isotope_element(53)="Dy"
    isotope_mass(53)=157.924412d0
    isotope_abundance(53)=0.1d0
    isotope_element(54)="Dy"
    isotope_mass(54)=159.925203d0
    isotope_abundance(54)=2.34d0
    isotope_element(55)="Dy"
    isotope_mass(55)=160.926939d0
    isotope_abundance(55)=18.9d0
    isotope_element(56)="Dy"
    isotope_mass(56)=161.926805d0
    isotope_abundance(56)=25.5d0
    isotope_element(57)="Dy"
    isotope_mass(57)=162.928737d0
    isotope_abundance(57)=24.9d0
    isotope_element(58)="Dy"
    isotope_mass(58)=163.929183d0
    isotope_abundance(58)=28.2d0
    isotope_element(59)="Er"
    isotope_mass(59)=161.928787d0
    isotope_abundance(59)=0.14d0
    isotope_element(60)="Er"
    isotope_mass(60)=163.929211d0
    isotope_abundance(60)=1.61d0
    isotope_element(61)="Er"
    isotope_mass(61)=165.930305d0
    isotope_abundance(61)=33.6d0
    isotope_element(62)="Er"
    isotope_mass(62)=166.932061d0
    isotope_abundance(62)=22.95d0
    isotope_element(63)="Er"
    isotope_mass(63)=167.932383d0
    isotope_abundance(63)=26.8d0
    isotope_element(64)="Er"
    isotope_mass(64)=169.935476d0
    isotope_abundance(64)=14.9d0
    isotope_element(65)="Eu"
    isotope_mass(65)=150.91986d0
    isotope_abundance(65)=47.8d0
    isotope_element(66)="Eu"
    isotope_mass(66)=152.921243d0
    isotope_abundance(66)=52.2d0
    isotope_element(67)="F"
    isotope_mass(67)=18.998403d0
    isotope_abundance(67)=100.0d0
    isotope_element(68)="Fe"
    isotope_mass(68)=53.939612d0
    isotope_abundance(68)=5.8d0
    isotope_element(69)="Fe"
    isotope_mass(69)=55.934939d0
    isotope_abundance(69)=91.72d0
    isotope_element(70)="Fe"
    isotope_mass(70)=56.935396d0
    isotope_abundance(70)=2.2d0
    isotope_element(71)="Fe"
    isotope_mass(71)=57.933278d0
    isotope_abundance(71)=0.28d0
    isotope_element(72)="Ga"
    isotope_mass(72)=68.925581d0
    isotope_abundance(72)=60.1d0
    isotope_element(73)="Ga"
    isotope_mass(73)=70.924701d0
    isotope_abundance(73)=39.9d0
    isotope_element(74)="Gd"
    isotope_mass(74)=151.919803d0
    isotope_abundance(74)=0.2d0
    isotope_element(75)="Gd"
    isotope_mass(75)=153.920876d0
    isotope_abundance(75)=2.18d0
    isotope_element(76)="Gd"
    isotope_mass(76)=154.822629d0
    isotope_abundance(76)=14.8d0
    isotope_element(77)="Gd"
    isotope_mass(77)=155.92213d0
    isotope_abundance(77)=20.47d0
    isotope_element(78)="Gd"
    isotope_mass(78)=156.923967d0
    isotope_abundance(78)=15.65d0
    isotope_element(79)="Gd"
    isotope_mass(79)=157.924111d0
    isotope_abundance(79)=24.84d0
    isotope_element(80)="Gd"
    isotope_mass(80)=159.927061d0
    isotope_abundance(80)=21.86d0
    isotope_element(81)="Ge"
    isotope_mass(81)=69.92425d0
    isotope_abundance(81)=20.5d0
    isotope_element(82)="Ge"
    isotope_mass(82)=71.92208d0
    isotope_abundance(82)=27.4d0
    isotope_element(83)="Ge"
    isotope_mass(83)=72.923464d0
    isotope_abundance(83)=7.8d0
    isotope_element(84)="Ge"
    isotope_mass(84)=73.921179d0
    isotope_abundance(84)=36.5d0
    isotope_element(85)="Ge"
    isotope_mass(85)=75.921403d0
    isotope_abundance(85)=7.8d0
    isotope_element(86)="H"
    isotope_mass(86)=1.007825d0
    isotope_abundance(86)=99.99d0
    isotope_element(87)="H"
    isotope_mass(87)=2.014102d0
    isotope_abundance(87)=0.015d0
    isotope_element(88)="He"
    isotope_mass(88)=3.016029d0
    isotope_abundance(88)=0.0001d0
    isotope_element(89)="He"
    isotope_mass(89)=4.002603d0
    isotope_abundance(89)=100.0d0
    isotope_element(90)="Hf"
    isotope_mass(90)=173.940065d0
    isotope_abundance(90)=0.16d0
    isotope_element(91)="Hf"
    isotope_mass(91)=175.94142d0
    isotope_abundance(91)=5.2d0
    isotope_element(92)="Hf"
    isotope_mass(92)=176.943233d0
    isotope_abundance(92)=18.6d0
    isotope_element(93)="Hf"
    isotope_mass(93)=177.94371d0
    isotope_abundance(93)=27.1d0
    isotope_element(94)="Hf"
    isotope_mass(94)=178.945827d0
    isotope_abundance(94)=13.74d0
    isotope_element(95)="Hf"
    isotope_mass(95)=179.946561d0
    isotope_abundance(95)=35.2d0
    isotope_element(96)="Hg"
    isotope_mass(96)=195.965812d0
    isotope_abundance(96)=0.15d0
    isotope_element(97)="Hg"
    isotope_mass(97)=197.96676d0
    isotope_abundance(97)=10.1d0
    isotope_element(98)="Hg"
    isotope_mass(98)=198.968269d0
    isotope_abundance(98)=17.0d0
    isotope_element(99)="Hg"
    isotope_mass(99)=199.968316d0
    isotope_abundance(99)=23.1d0
    isotope_element(100)="Hg"
    isotope_mass(100)=200.970293d0
    isotope_abundance(100)=13.2d0
    isotope_element(101)="Hg"
    isotope_mass(101)=201.970632d0
    isotope_abundance(101)=29.65d0
    isotope_element(102)="Hg"
    isotope_mass(102)=203.973481d0
    isotope_abundance(102)=6.8d0
    isotope_element(103)="Ho"
    isotope_mass(103)=164.930332d0
    isotope_abundance(103)=100.0d0
    isotope_element(104)="I"
    isotope_mass(104)=126.904477d0
    isotope_abundance(104)=100.0d0
    isotope_element(105)="In"
    isotope_mass(105)=112.904056d0
    isotope_abundance(105)=4.3d0
    isotope_element(106)="In"
    isotope_mass(106)=114.903875d0
    isotope_abundance(106)=95.7d0
    isotope_element(107)="Ir"
    isotope_mass(107)=190.960603d0
    isotope_abundance(107)=37.3d0
    isotope_element(108)="Ir"
    isotope_mass(108)=192.962942d0
    isotope_abundance(108)=62.7d0
    isotope_element(109)="K"
    isotope_mass(109)=38.963708d0
    isotope_abundance(109)=93.2d0
    isotope_element(110)="K"
    isotope_mass(110)=39.963999d0
    isotope_abundance(110)=0.012d0
    isotope_element(111)="K"
    isotope_mass(111)=40.961825d0
    isotope_abundance(111)=6.73d0
    isotope_element(112)="Kr"
    isotope_mass(112)=77.920397d0
    isotope_abundance(112)=0.35d0
    isotope_element(113)="Kr"
    isotope_mass(113)=79.916375d0
    isotope_abundance(113)=2.25d0
    isotope_element(114)="Kr"
    isotope_mass(114)=81.913483d0
    isotope_abundance(114)=11.6d0
    isotope_element(115)="Kr"
    isotope_mass(115)=82.914134d0
    isotope_abundance(115)=11.5d0
    isotope_element(116)="Kr"
    isotope_mass(116)=83.911506d0
    isotope_abundance(116)=57.0d0
    isotope_element(117)="Kr"
    isotope_mass(117)=85.910614d0
    isotope_abundance(117)=17.3d0
    isotope_element(118)="La"
    isotope_mass(118)=137.907114d0
    isotope_abundance(118)=0.09d0
    isotope_element(119)="La"
    isotope_mass(119)=138.906355d0
    isotope_abundance(119)=99.91d0
    isotope_element(120)="Li"
    isotope_mass(120)=6.015123d0
    isotope_abundance(120)=7.42d0
    isotope_element(121)="Li"
    isotope_mass(121)=7.016005d0
    isotope_abundance(121)=92.58d0
    isotope_element(122)="Lu"
    isotope_mass(122)=174.940785d0
    isotope_abundance(122)=97.4d0
    isotope_element(123)="Lu"
    isotope_mass(123)=175.942694d0
    isotope_abundance(123)=2.6d0
    isotope_element(124)="Mg"
    isotope_mass(124)=23.985045d0
    isotope_abundance(124)=78.9d0
    isotope_element(125)="Mg"
    isotope_mass(125)=24.985839d0
    isotope_abundance(125)=10.0d0
    isotope_element(126)="Mg"
    isotope_mass(126)=25.982595d0
    isotope_abundance(126)=11.1d0
    isotope_element(127)="Mn"
    isotope_mass(127)=54.938046d0
    isotope_abundance(127)=100.0d0
    isotope_element(128)="Mo"
    isotope_mass(128)=91.906809d0
    isotope_abundance(128)=14.84d0
    isotope_element(129)="Mo"
    isotope_mass(129)=93.905086d0
    isotope_abundance(129)=9.25d0
    isotope_element(130)="Mo"
    isotope_mass(130)=94.905838d0
    isotope_abundance(130)=15.92d0
    isotope_element(131)="Mo"
    isotope_mass(131)=95.904676d0
    isotope_abundance(131)=16.68d0
    isotope_element(132)="Mo"
    isotope_mass(132)=96.906018d0
    isotope_abundance(132)=9.55d0
    isotope_element(133)="Mo"
    isotope_mass(133)=97.905405d0
    isotope_abundance(133)=24.13d0
    isotope_element(134)="Mo"
    isotope_mass(134)=99.907473d0
    isotope_abundance(134)=9.63d0
    isotope_element(135)="N"
    isotope_mass(135)=14.003074d0
    isotope_abundance(135)=99.63d0
    isotope_element(136)="N"
    isotope_mass(136)=15.000109d0
    isotope_abundance(136)=0.37d0
    isotope_element(137)="Na"
    isotope_mass(137)=22.98977d0
    isotope_abundance(137)=100.0d0
    isotope_element(138)="Nb"
    isotope_mass(138)=92.906378d0
    isotope_abundance(138)=100.0d0
    isotope_element(139)="Nd"
    isotope_mass(139)=141.907731d0
    isotope_abundance(139)=27.13d0
    isotope_element(140)="Nd"
    isotope_mass(140)=142.909823d0
    isotope_abundance(140)=12.18d0
    isotope_element(141)="Nd"
    isotope_mass(141)=143.910096d0
    isotope_abundance(141)=23.8d0
    isotope_element(142)="Nd"
    isotope_mass(142)=144.912582d0
    isotope_abundance(142)=8.3d0
    isotope_element(143)="Nd"
    isotope_mass(143)=145.913126d0
    isotope_abundance(143)=17.19d0
    isotope_element(144)="Nd"
    isotope_mass(144)=147.916901d0
    isotope_abundance(144)=5.76d0
    isotope_element(145)="Nd"
    isotope_mass(145)=149.9209d0
    isotope_abundance(145)=5.64d0
    isotope_element(146)="Ne"
    isotope_mass(146)=19.992439d0
    isotope_abundance(146)=90.6d0
    isotope_element(147)="Ne"
    isotope_mass(147)=20.993845d0
    isotope_abundance(147)=0.26d0
    isotope_element(148)="Ne"
    isotope_mass(148)=21.991384d0
    isotope_abundance(148)=9.2d0
    isotope_element(149)="Ni"
    isotope_mass(149)=57.935347d0
    isotope_abundance(149)=68.27d0
    isotope_element(150)="Ni"
    isotope_mass(150)=59.930789d0
    isotope_abundance(150)=26.1d0
    isotope_element(151)="Ni"
    isotope_mass(151)=60.931059d0
    isotope_abundance(151)=1.13d0
    isotope_element(152)="Ni"
    isotope_mass(152)=61.928346d0
    isotope_abundance(152)=3.59d0
    isotope_element(153)="Ni"
    isotope_mass(153)=63.927968d0
    isotope_abundance(153)=0.91d0
    isotope_element(154)="O"
    isotope_mass(154)=15.994915d0
    isotope_abundance(154)=99.76d0
    isotope_element(155)="O"
    isotope_mass(155)=16.999131d0
    isotope_abundance(155)=0.038d0
    isotope_element(156)="O"
    isotope_mass(156)=17.999159d0
    isotope_abundance(156)=0.2d0
    isotope_element(157)="Os"
    isotope_mass(157)=183.952514d0
    isotope_abundance(157)=0.02d0
    isotope_element(158)="Os"
    isotope_mass(158)=185.953852d0
    isotope_abundance(158)=1.58d0
    isotope_element(159)="Os"
    isotope_mass(159)=186.955762d0
    isotope_abundance(159)=1.6d0
    isotope_element(160)="Os"
    isotope_mass(160)=187.95585d0
    isotope_abundance(160)=13.3d0
    isotope_element(161)="Os"
    isotope_mass(161)=188.958156d0
    isotope_abundance(161)=16.1d0
    isotope_element(162)="Os"
    isotope_mass(162)=189.958455d0
    isotope_abundance(162)=26.4d0
    isotope_element(163)="Os"
    isotope_mass(163)=191.961487d0
    isotope_abundance(163)=41.0d0
    isotope_element(164)="P"
    isotope_mass(164)=30.973763d0
    isotope_abundance(164)=100.0d0
    isotope_element(165)="Pb"
    isotope_mass(165)=203.973037d0
    isotope_abundance(165)=1.4d0
    isotope_element(166)="Pb"
    isotope_mass(166)=205.974455d0
    isotope_abundance(166)=24.1d0
    isotope_element(167)="Pb"
    isotope_mass(167)=206.975885d0
    isotope_abundance(167)=22.1d0
    isotope_element(168)="Pb"
    isotope_mass(168)=207.976641d0
    isotope_abundance(168)=52.4d0
    isotope_element(169)="Pd"
    isotope_mass(169)=101.905609d0
    isotope_abundance(169)=1.02d0
    isotope_element(170)="Pd"
    isotope_mass(170)=103.904026d0
    isotope_abundance(170)=11.14d0
    isotope_element(171)="Pd"
    isotope_mass(171)=104.905075d0
    isotope_abundance(171)=22.33d0
    isotope_element(172)="Pd"
    isotope_mass(172)=105.903475d0
    isotope_abundance(172)=27.33d0
    isotope_element(173)="Pd"
    isotope_mass(173)=107.903894d0
    isotope_abundance(173)=26.46d0
    isotope_element(174)="Pd"
    isotope_mass(174)=109.905169d0
    isotope_abundance(174)=11.72d0
    isotope_element(175)="Pr"
    isotope_mass(175)=140.907657d0
    isotope_abundance(175)=100.0d0
    isotope_element(176)="Pt"
    isotope_mass(176)=189.959937d0
    isotope_abundance(176)=0.01d0
    isotope_element(177)="Pt"
    isotope_mass(177)=191.961049d0
    isotope_abundance(177)=0.79d0
    isotope_element(178)="Pt"
    isotope_mass(178)=193.962679d0
    isotope_abundance(178)=32.9d0
    isotope_element(179)="Pt"
    isotope_mass(179)=194.964785d0
    isotope_abundance(179)=33.8d0
    isotope_element(180)="Pt"
    isotope_mass(180)=195.964947d0
    isotope_abundance(180)=25.3d0
    isotope_element(181)="Pt"
    isotope_mass(181)=197.967879d0
    isotope_abundance(181)=7.2d0
    isotope_element(182)="Rb"
    isotope_mass(182)=84.9118d0
    isotope_abundance(182)=72.17d0
    isotope_element(183)="Rb"
    isotope_mass(183)=86.909184d0
    isotope_abundance(183)=27.84d0
    isotope_element(184)="Re"
    isotope_mass(184)=184.952977d0
    isotope_abundance(184)=37.4d0
    isotope_element(185)="Re"
    isotope_mass(185)=186.955765d0
    isotope_abundance(185)=62.6d0
    isotope_element(186)="Rh"
    isotope_mass(186)=102.905503d0
    isotope_abundance(186)=100.0d0
    isotope_element(187)="Ru"
    isotope_mass(187)=95.907596d0
    isotope_abundance(187)=5.52d0
    isotope_element(188)="Ru"
    isotope_mass(188)=97.905287d0
    isotope_abundance(188)=1.88d0
    isotope_element(189)="Ru"
    isotope_mass(189)=98.905937d0
    isotope_abundance(189)=12.7d0
    isotope_element(190)="Ru"
    isotope_mass(190)=99.904218d0
    isotope_abundance(190)=12.6d0
    isotope_element(191)="Ru"
    isotope_mass(191)=100.905581d0
    isotope_abundance(191)=17.0d0
    isotope_element(192)="Ru"
    isotope_mass(192)=101.904348d0
    isotope_abundance(192)=31.6d0
    isotope_element(193)="Ru"
    isotope_mass(193)=103.905422d0
    isotope_abundance(193)=18.7d0
    isotope_element(194)="S"
    isotope_mass(194)=31.972072d0
    isotope_abundance(194)=95.02d0
    isotope_element(195)="S"
    isotope_mass(195)=32.971459d0
    isotope_abundance(195)=0.75d0
    isotope_element(196)="S"
    isotope_mass(196)=33.967868d0
    isotope_abundance(196)=4.21d0
    isotope_element(197)="S"
    isotope_mass(197)=35.967079d0
    isotope_abundance(197)=0.02d0
    isotope_element(198)="Sb"
    isotope_mass(198)=120.903824d0
    isotope_abundance(198)=57.3d0
    isotope_element(199)="Sb"
    isotope_mass(199)=122.904222d0
    isotope_abundance(199)=42.7d0
    isotope_element(200)="Sc"
    isotope_mass(200)=44.955914d0
    isotope_abundance(200)=100.0d0
    isotope_element(201)="Se"
    isotope_mass(201)=73.922477d0
    isotope_abundance(201)=0.9d0
    isotope_element(202)="Se"
    isotope_mass(202)=75.919207d0
    isotope_abundance(202)=9.0d0
    isotope_element(203)="Se"
    isotope_mass(203)=76.919908d0
    isotope_abundance(203)=7.6d0
    isotope_element(204)="Se"
    isotope_mass(204)=77.917304d0
    isotope_abundance(204)=23.5d0
    isotope_element(205)="Se"
    isotope_mass(205)=79.916521d0
    isotope_abundance(205)=49.6d0
    isotope_element(206)="Se"
    isotope_mass(206)=81.916709d0
    isotope_abundance(206)=9.4d0
    isotope_element(207)="Si"
    isotope_mass(207)=27.976928d0
    isotope_abundance(207)=92.23d0
    isotope_element(208)="Si"
    isotope_mass(208)=28.976496d0
    isotope_abundance(208)=4.67d0
    isotope_element(209)="Si"
    isotope_mass(209)=29.973772d0
    isotope_abundance(209)=3.1d0
    isotope_element(210)="Sm"
    isotope_mass(210)=143.912009d0
    isotope_abundance(210)=3.1d0
    isotope_element(211)="Sm"
    isotope_mass(211)=146.914907d0
    isotope_abundance(211)=15.0d0
    isotope_element(212)="Sm"
    isotope_mass(212)=147.914832d0
    isotope_abundance(212)=11.3d0
    isotope_element(213)="Sm"
    isotope_mass(213)=148.917193d0
    isotope_abundance(213)=13.8d0
    isotope_element(214)="Sm"
    isotope_mass(214)=149.917285d0
    isotope_abundance(214)=7.4d0
    isotope_element(215)="Sm"
    isotope_mass(215)=151.919741d0
    isotope_abundance(215)=26.7d0
    isotope_element(216)="Sm"
    isotope_mass(216)=153.922218d0
    isotope_abundance(216)=22.7d0
    isotope_element(217)="Sn"
    isotope_mass(217)=111.904826d0
    isotope_abundance(217)=0.97d0
    isotope_element(218)="Sn"
    isotope_mass(218)=113.902784d0
    isotope_abundance(218)=0.65d0
    isotope_element(219)="Sn"
    isotope_mass(219)=114.903348d0
    isotope_abundance(219)=0.36d0
    isotope_element(220)="Sn"
    isotope_mass(220)=115.901744d0
    isotope_abundance(220)=14.7d0
    isotope_element(221)="Sn"
    isotope_mass(221)=116.902954d0
    isotope_abundance(221)=7.7d0
    isotope_element(222)="Sn"
    isotope_mass(222)=117.901607d0
    isotope_abundance(222)=24.3d0
    isotope_element(223)="Sn"
    isotope_mass(223)=118.90331d0
    isotope_abundance(223)=8.6d0
    isotope_element(224)="Sn"
    isotope_mass(224)=119.902199d0
    isotope_abundance(224)=32.4d0
    isotope_element(225)="Sn"
    isotope_mass(225)=121.90344d0
    isotope_abundance(225)=4.6d0
    isotope_element(226)="Sn"
    isotope_mass(226)=123.905271d0
    isotope_abundance(226)=5.6d0
    isotope_element(227)="Sr"
    isotope_mass(227)=83.913428d0
    isotope_abundance(227)=0.56d0
    isotope_element(228)="Sr"
    isotope_mass(228)=85.909273d0
    isotope_abundance(228)=9.86d0
    isotope_element(229)="Sr"
    isotope_mass(229)=86.908902d0
    isotope_abundance(229)=7.0d0
    isotope_element(230)="Sr"
    isotope_mass(230)=87.905625d0
    isotope_abundance(230)=82.58d0
    isotope_element(231)="Ta"
    isotope_mass(231)=179.947489d0
    isotope_abundance(231)=0.012d0
    isotope_element(232)="Ta"
    isotope_mass(232)=180.948014d0
    isotope_abundance(232)=99.99d0
    isotope_element(233)="Tb"
    isotope_mass(233)=158.92535d0
    isotope_abundance(233)=100.0d0
    isotope_element(234)="Te"
    isotope_mass(234)=119.904021d0
    isotope_abundance(234)=0.096d0
    isotope_element(235)="Te"
    isotope_mass(235)=121.903055d0
    isotope_abundance(235)=2.6d0
    isotope_element(236)="Te"
    isotope_mass(236)=122.904278d0
    isotope_abundance(236)=0.91d0
    isotope_element(237)="Te"
    isotope_mass(237)=123.902825d0
    isotope_abundance(237)=4.82d0
    isotope_element(238)="Te"
    isotope_mass(238)=124.904435d0
    isotope_abundance(238)=7.14d0
    isotope_element(239)="Te"
    isotope_mass(239)=125.90331d0
    isotope_abundance(239)=18.95d0
    isotope_element(240)="Te"
    isotope_mass(240)=127.904464d0
    isotope_abundance(240)=31.69d0
    isotope_element(241)="Te"
    isotope_mass(241)=129.906229d0
    isotope_abundance(241)=33.8d0
    isotope_element(242)="Th"
    isotope_mass(242)=232.038054d0
    isotope_abundance(242)=100.0d0
    isotope_element(243)="Ti"
    isotope_mass(243)=45.952633d0
    isotope_abundance(243)=8.0d0
    isotope_element(244)="Ti"
    isotope_mass(244)=46.951765d0
    isotope_abundance(244)=7.3d0
    isotope_element(245)="Ti"
    isotope_mass(245)=47.947947d0
    isotope_abundance(245)=73.8d0
    isotope_element(246)="Ti"
    isotope_mass(246)=48.947871d0
    isotope_abundance(246)=5.5d0
    isotope_element(247)="Ti"
    isotope_mass(247)=49.944786d0
    isotope_abundance(247)=5.4d0
    isotope_element(248)="Tl"
    isotope_mass(248)=202.972336d0
    isotope_abundance(248)=29.52d0
    isotope_element(249)="Tl"
    isotope_mass(249)=204.97441d0
    isotope_abundance(249)=70.48d0
    isotope_element(250)="Tm"
    isotope_mass(250)=168.934225d0
    isotope_abundance(250)=100.0d0
    isotope_element(251)="U"
    isotope_mass(251)=234.040947d0
    isotope_abundance(251)=0.006d0
    isotope_element(252)="U"
    isotope_mass(252)=235.043925d0
    isotope_abundance(252)=0.72d0
    isotope_element(253)="U"
    isotope_mass(253)=238.050786d0
    isotope_abundance(253)=99.27d0
    isotope_element(254)="V"
    isotope_mass(254)=49.947161d0
    isotope_abundance(254)=0.25d0
    isotope_element(255)="V"
    isotope_mass(255)=50.943963d0
    isotope_abundance(255)=99.75d0
    isotope_element(256)="W"
    isotope_mass(256)=179.946727d0
    isotope_abundance(256)=0.13d0
    isotope_element(257)="W"
    isotope_mass(257)=181.948225d0
    isotope_abundance(257)=26.3d0
    isotope_element(258)="W"
    isotope_mass(258)=182.950245d0
    isotope_abundance(258)=14.3d0
    isotope_element(259)="W"
    isotope_mass(259)=183.950953d0
    isotope_abundance(259)=30.67d0
    isotope_element(260)="W"
    isotope_mass(260)=185.954377d0
    isotope_abundance(260)=28.6d0
    isotope_element(261)="Xe"
    isotope_mass(261)=123.905894d0
    isotope_abundance(261)=0.1d0
    isotope_element(262)="Xe"
    isotope_mass(262)=125.904281d0
    isotope_abundance(262)=0.09d0
    isotope_element(263)="Xe"
    isotope_mass(263)=127.903531d0
    isotope_abundance(263)=1.91d0
    isotope_element(264)="Xe"
    isotope_mass(264)=128.90478d0
    isotope_abundance(264)=26.4d0
    isotope_element(265)="Xe"
    isotope_mass(265)=129.90351d0
    isotope_abundance(265)=4.1d0
    isotope_element(266)="Xe"
    isotope_mass(266)=130.905076d0
    isotope_abundance(266)=21.2d0
    isotope_element(267)="Xe"
    isotope_mass(267)=131.904148d0
    isotope_abundance(267)=26.9d0
    isotope_element(268)="Xe"
    isotope_mass(268)=133.905395d0
    isotope_abundance(268)=10.4d0
    isotope_element(269)="Xe"
    isotope_mass(269)=135.907219d0
    isotope_abundance(269)=8.9d0
    isotope_element(270)="Y"
    isotope_mass(270)=88.905856d0
    isotope_abundance(270)=100.0d0
    isotope_element(271)="Yb"
    isotope_mass(271)=167.933908d0
    isotope_abundance(271)=0.13d0
    isotope_element(272)="Yb"
    isotope_mass(272)=169.934774d0
    isotope_abundance(272)=3.05d0
    isotope_element(273)="Yb"
    isotope_mass(273)=170.936338d0
    isotope_abundance(273)=14.3d0
    isotope_element(274)="Yb"
    isotope_mass(274)=171.936393d0
    isotope_abundance(274)=21.9d0
    isotope_element(275)="Yb"
    isotope_mass(275)=172.938222d0
    isotope_abundance(275)=16.12d0
    isotope_element(276)="Yb"
    isotope_mass(276)=173.938873d0
    isotope_abundance(276)=31.8d0
    isotope_element(277)="Yb"
    isotope_mass(277)=175.942576d0
    isotope_abundance(277)=12.7d0
    isotope_element(278)="Zn"
    isotope_mass(278)=63.929145d0
    isotope_abundance(278)=48.6d0
    isotope_element(279)="Zn"
    isotope_mass(279)=65.926035d0
    isotope_abundance(279)=27.9d0
    isotope_element(280)="Zn"
    isotope_mass(280)=66.927129d0
    isotope_abundance(280)=4.1d0
    isotope_element(281)="Zn"
    isotope_mass(281)=67.924846d0
    isotope_abundance(281)=18.8d0
    isotope_element(282)="Zn"
    isotope_mass(282)=69.925325d0
    isotope_abundance(282)=0.6d0
    isotope_element(283)="Zr"
    isotope_mass(283)=89.904708d0
    isotope_abundance(283)=51.45d0
    isotope_element(284)="Zr"
    isotope_mass(284)=90.905644d0
    isotope_abundance(284)=11.27d0
    isotope_element(285)="Zr"
    isotope_mass(285)=91.905039d0
    isotope_abundance(285)=17.17d0
    isotope_element(286)="Zr"
    isotope_mass(286)=93.906319d0
    isotope_abundance(286)=17.33d0
    isotope_element(287)="Zr"
    isotope_mass(287)=95.908272d0
    isotope_abundance(287)=2.78d0
  end subroutine data_fill_isotopes

  ! Free the memory used by the isotopic information.
  subroutine data_free_isotopes()
    deallocate(isotope_element)
    deallocate(isotope_mass)
    deallocate(isotope_abundance)
  end subroutine data_free_isotopes

  ! Compute the average mass of each element and its g-factor (Pearson
  ! deviation coefficient of the masses).
  subroutine data_calc_mandg(element,m,g)
    character(len=3),intent(in) :: element
    real(kind=8),intent(out) :: m,g

    integer(kind=4) :: i,niso

    niso=size(isotope_element)

    m=0.
    do i=1,niso
       if(isotope_element(i).eq.element) then
          m=m+isotope_mass(i)*isotope_abundance(i)
       end if
    end do
    m=m/100.
    g=0.
    do i=1,niso
       if(isotope_element(i).eq.element) then
          g=g+isotope_abundance(i)*(1.-isotope_mass(i)/m)**2
       end if
    end do
    g=g/100.
  end subroutine data_calc_mandg
end module data
