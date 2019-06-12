module shared_arr
double precision, allocatable::t_arr(:),t(:,:),w(:),z(:),ep(:),test(:)
complex*16 ,allocatable::up_up(:,:),up_down(:,:),down_up(:,:),down_down(:,:)
double precision ,allocatable::n_a_up(:),n_a_dwn(:),n_b_up(:),n_b_dwn(:),num(:)
double precision, allocatable::bincenergy(:),rwork(:)
double precision, allocatable::bincnumber(:),prob(:),np(:)
complex*16, allocatable::phi_x(:),zhi_x(:),phi_y(:),zhi_y(:),phi_z(:),zhi_z(:)
double precision, allocatable::dis(:)
complex*16, allocatable::H(:,:),work(:)
integer, allocatable :: pi(:),pj(:)
complex*16, allocatable :: Sz_local(:,:), Sx_local(:,:), Sy_local(:,:)
complex*16, allocatable :: Tz_local(:,:), Tx_local(:,:), Ty_local(:,:)
complex*16, allocatable :: SP_Density_Matrix(:,:),S_i_S_j(:,:),T_i_T_j(:,:)
end module

module glob_var
	integer shift,m_size,fl,iter,count,flag_t,flag_mu,d,seed,seedC,seedV,orbital
	double precision temp,mus,np1,alpha
	double precision n_av,n_av2,muav,mu_sum2,av_count,r
    double precision :: U,J_val,m_d,fill,V,filling
    double precision :: A1,B1,e_onsite 
    complex*16       :: AA1,AA2,AA3,AA4,AA5,AA6,AA7,AA8,AA9,AA10,AA11,AA12,i_complex
    double precision :: t_aa,t_ab,t_ba,t_bb,ipr2,ipr_a,ipr_b
    complex*16 :: osf,ssf,qosf,qsf
    complex*16 :: qof_xy1, qof_xy2, qof_xy3, qof_xy4, qof_xy5, qof_xy6, qof_xy7, qof_xy8
    complex*16 :: qof_z1, qof_z2, qof_z3, qof_z4
    complex*16 :: qsf_xy1, qsf_xy2, qsf_xy3, qsf_xy4
    complex*16 :: qsf_z1, qsf_z2, qsf_z3, qsf_z4, qsf_z5, qsf_z6, qsf_z7, qsf_z8, qsf_z9, qsf_z10, qsf_z11, qsf_z12
    complex*16 :: Quant_orb_strfac,Quant_spin_strfac
end module


!************************ Main program starts*************************************************

program hubbard
use mtmod
use shared_arr
use glob_var
double precision :: orb,binsize,binsize_n,gama,pii
double precision :: f0, f, fL2, fR, m, mR, mL, rtmp,f_t
double precision :: energy1, energy2,energy,energy_cl1, energy_cl2,chem_pot
double precision :: BB1,q,numbr,dos_en,dos_en1,gama_d,delta_mu
double precision ::BB2,input
double precision :: op1,op2,op3,op4,op5,op6
complex*16 ::dos_aa_up,dos_bb_up,dos_aa_dwn,dos_bb_dwn,dos
double precision :: BB3,BB4,BB5,BB6,q_x,q_y
double precision ::en_mu_l,en_mu_u,prob_sum
integer ::i,j,k,p,ww,l,omega,kk,i_nx,i_ny
integer :: id,ii,jd,ji,i_i,j_j,a,b,nbin,nn
integer ::iii,count1,ref,sum1,count2,site
integer :: int_read_seed


    open(unit=12,file="input.inp",status="unknown")

    do i=1,18
    read(12,*) input
    if(i.eq.1)orbital=int(input)
    if(i.eq.2)d=int(input)
    if(i.eq.3)alpha= dble(input)
    if(i.eq.4)U=dble(input)
    if(i.eq.5)filling =dble(input)
    if(i.eq.6)temp=dble(input)
    if(i.eq.7)seedC=int(input)
    if(i.eq.8)seedV=int(input)
    if(i.eq.9) V=dble(input)
    if(i.eq.10)t_aa=dble(input)
    if(i.eq.11)t_ab=dble(input)
    if(i.eq.12)t_ba=dble(input)
    if(i.eq.13)t_bb=dble(input)
    if(i.eq.14)delta_mu=dble(input)
    if(i.eq.15)gama=dble(input)
    if(i.eq.16)gama_d=dble(input)
    if(i.eq.17)J_val=dble(input)
    if(i.eq.18)int_read_seed=int(input) 
    enddo
    close(12)
     
    print*,'orbital=',orbital
    print*,'dimension=',d
    print*,'alpha=',alpha
    print*,'U=',U
    print*,'filling=',filling
    print*,'temp=',temp
    print*,'config seed=',seedC
    print*,'disorder seed=',seedV
    print*,'disorder=',V
    print*,'t_aa=',t_aa
    print*,'t_ab=',t_ab
    print*,'t_ba=',t_ba
    print*,'t_bb=',t_bb
    print*,'delta_mu',delta_mu
    print*,'gama_p(n)',gama  
    print*,'gama_dos',gama_d 
    print*,'J',J_val
    print*,'int_read_seed', int_read_seed



shift = orbital
pii = acos(-1.0)
!J_val = 0.0d0*U/4.0d0
m_size =4*(d**2)            !Matrix size!
fill = 2*d**2*filling       !Filling (no. of particles)!
count1 = 0
A1 = 1.0d0*(U-2*J_val)
B1 = 1.0d0*(U-3*J_val)
print*, "A1=",A1,"   B1=",B1
i_complex = cmplx(0.0d0,1.0d0)
allocate (t_arr(shift**2), t(shift,shift),w(m_size),z(nmax))
allocate (up_up(d**2,shift**2),up_down(d**2,shift**2),test(m_size))
allocate (down_up(d**2,shift**2),down_down(d**2,shift**2))
allocate (H(m_size,m_size))
allocate (bincenergy(nbin),ep(m_size),pi(d**2),pj(d**2))
allocate (bincnumber(nbin),np(d**2),prob(500))
allocate (work(2*4*d**2-1),rwork(3*4*d**2-2))
allocate (phi_x(d**2),zhi_x(d**2),phi_y(d**2),zhi_y(d**2),phi_z(d**2),zhi_z(d**2),dis(d**2))
allocate (n_a_up(d**2),n_a_dwn(d**2),n_b_up(d**2),n_b_dwn(d**2),num(d**2))
allocate (SP_Density_Matrix(4*d**2,4*d**2))  
allocate (Sz_local(4,4),Sy_local(4,4),Sx_local(4,4))
allocate (Tz_local(4,4),Ty_local(4,4),Tx_local(4,4))
allocate (S_i_S_j(d**2,d**2),T_i_T_j(d**2,d**2))

open (14,file = 'MF_param_biased.dat')
open (15, file ='consis_steps.dat')
open (26, file = 'energy_file.dat')
open (27, file ='num_densities.dat')
open (28, file ='Output_OrderParams.dat')
open (33, file ='eigen_val.dat')
open (31, file ='ipr.dat')
open (39, file ='dis_prof.dat')
open (41, file = 'ipr_a.dat')
open (42, file = 'ipr_b.dat')
open (61, file ='ipr_c.dat')
open (71, file = 'ipr_ac.dat')
open (72, file = 'ipr_bc.dat')
open (85, file = 'spin_sf.dat')
open (86, file = 'orb_sf.dat')
open (90, file = 'qspin_sf.dat')
open (95, file = 'qorb_sf.dat')
open (100,file ='dos.dat')
open (101,file ='prob_n.dat')
open (102,file ='dos_a.dat')
open (103,file ='dos_b.dat')
open (104,file ='dos_a+b.dat')
!open (185,file ='qosf.dat')
open (200,file ='local_S2_T2.dat')
open (1001, file ='Initial_OrderParams.dat')
open (1002 ,file = 'Input_OrderParams.dat')



!Matrix label generator

 call sgrnd(seedV)
do i =1,d**2
   call rannum(r)
   dis(i) = 2.0d0*(r-0.50d0)*V
write(39,*) i,dis(i)
 enddo


     

do l=1,d**2
          do p=1,d
           do ww=1,d
            if (d*(p-1)+ww.eq.l) then
             pi(l)=p
             pj(l)=ww
             go to 50
            end if
           end do
50       end do
end do


!************************ Initial conditions for phi_z() and zhi_z() *************************
!Condition 1 : Staggered config
goto 11
 do i =1,d**2
    if (mod((pi(i)+pj(i)),2).eq. 0) then
       phi_x(i) = 1.0d0
       zhi_x(i) =1.0d0
       phi_y(i) = 1.0d0
       zhi_y(i) =1.0d0
       phi_z(i) = 1.0d0
       zhi_z(i) =1.0d0
    else
       phi_x(i) = -1.0d0
       zhi_x(i) =-1.0d0
       phi_y(i) = -1.0d0
       zhi_y(i) =-1.0d0
      phi_z(i) =-1.0d0
      zhi_z(i) = -1.0d0
    endif
 enddo
11 continue
!Condition 2: Random config

!goto 12
 
  call sgrnd(seedC)

if (int_read_seed .ne. 1) then
 do i =1,d**2
 call rannum(r)
 !write(*,*) r,'r'
 if (r.le. 0.50d0) then
 
 call rannum(r)

 phi_x(i) = r
 call rannum(r)
  zhi_x(i) = r
 call rannum(r)
  phi_y(i) = r
 call rannum(r)
  zhi_y(i) =r
 call rannum(r)   
 phi_z(i) = r
 call rannum(r)
 zhi_z(i) = r
 else
 call rannum(r)
 phi_x(i) = r
 call rannum(r)
 zhi_x(i) =r
 call rannum(r)
 phi_y(i) = r
 call rannum(r)
 zhi_y(i) =r
 call rannum(r)
 phi_z(i) = r
 call rannum(r)
 zhi_z(i) = r
 endif
 enddo
endif
!12 continue


 ! pause
 
!Condition 3: Initializing to some particular value

 !phi_z = 0.50d0
 !zhi_z = 0.50d0 



!Biased initial configuration
 
 if (int_read_seed .eq. 1) then
 do i =1,d**2
 read (1002,*) site,op1,op2,op3,op4,op5,op6
 phi_x(i) = cmplx(op1, 0.0)
 phi_y(i) = cmplx(op2, 0.0)
 phi_z(i) = cmplx(op3, 0.0)
 zhi_x(i) = cmplx(op4, 0.0)
 zhi_y(i) = cmplx(op5, 0.0)
 zhi_z(i) = cmplx(op6, 0.0)
 enddo
 endif


do i =1,d**2
 write (1001,*) i, real(phi_x(i)), real(phi_y(i)),  real(phi_z(i)), real(zhi_x(i)), real(zhi_y(i)), real(zhi_z(i))
 enddo
!*********************************************************************************************

 up_up=0.0d0
 up_down=0.0d0
 down_up=0.0d0
 down_down=0.0d0
        
        call matrix_gen
        call get_mu
        call mf

        energy = 0.0d0
        do i=1,m_size
     	  energy=energy+w(i)*(1.0d0/(1.0d0+exp((w(i)-m_d)/temp)))
	    enddo

      
       energy_cl1 = 0.0d0
      do i = 1, d**2
      energy_cl1 =  energy_cl1 + A1*(phi_x(i)**2+phi_y(i)**2+phi_z(i)**2)+&
                             B1*(zhi_x(i)**2+zhi_y(i)**2+zhi_z(i)**2)+&
                           ((3.0d0/4.0d0)*(2*U -3*J_val)*(n_a_up(i)+n_a_dwn(i)+n_b_up(i)+n_b_dwn(i)))
       enddo
        energy1 =energy+energy_cl1
       

       
       call matrix_gen
       call get_mu
       call mf
     
       energy2 = 0.0d0
       energy_cl2 = 0.0d0
       do i=1,m_size
       energy2=energy2+w(i)*(1.0d0/(1.0d0+exp((w(i)-m_d)/temp)))
	   enddo

       do i = 1, d**2
       energy_cl2 =  energy_cl2 + A1*(phi_x(i)**2+phi_y(i)**2+phi_z(i)**2)+&
                             B1*(zhi_x(i)**2+zhi_y(i)**2+zhi_z(i)**2)+&
                           ((3.0d0/4.0d0)*(2*U -3*J_val)*(n_a_up(i)+n_a_dwn(i)+n_b_up(i)+n_b_dwn(i)))
       enddo
        energy2 = energy2 + energy_cl2                  
	

	write(26,*) energy1,energy2, energy_cl1, energy_cl2,  (energy2 - energy1), 0



!***************************Consistency check*************************************************     
    
   
    do while ((abs(energy2-energy1).ge.1e-5).and. (count1.le.10000)) 
    
    
    energy1 = energy2
    energy_cl1=energy_cl2
    
    call mf
    call matrix_gen
    call get_mu

    
	
    energy2 = 0.0d0
    do i=1,m_size
       energy2=energy2+w(i)*(1/(1+exp((w(i)-m_d)/temp)))
    enddo
     energy_cl2 = 0.0d0
     do i = 1, d**2
      energy_cl2 =  energy_cl2 + A1*(phi_x(i)**2+phi_y(i)**2+phi_z(i)**2)+&
                             B1*(zhi_x(i)**2+zhi_y(i)**2+zhi_z(i)**2)+&
                           ((3.0d0/4.0d0)*(2*U -3*J_val)*(n_a_up(i)+n_a_dwn(i)+n_b_up(i)+n_b_dwn(i)))
       enddo
        energy2 = energy2 + energy_cl2                  
	
    count1 = count1+1
    write(26,*) energy1,energy2, energy_cl1, energy_cl2, (energy2 - energy1),count1
    call flush(26)
    f=0.0d0
	do i=1,m_size
	f=f+(1.0d0/(exp((w(i)-m_d)/temp)+1.0d0))
	end do

	write(14,*)f,m_d
    
    sum2 = 0.0d0
    !do i =1,d**2
    !write(10,*) i,(n_a_up(i)+n_a_dwn(i)+n_b_up(i)+n_b_dwn(i))
    !sum2 = sum2 + abs(n_a_up(i)+n_b_up(i)+n_a_dwn(i)+n_b_dwn(i))  !Total no.of particles!
    !enddo
    !write(15,*) count1,sum2
    !call flush(15)    


     enddo   

!********************************************************************************************
      !write(15,*) count1, 'count'
      
      call mf
      write(27,*) '#site','-------','n_a_up','-------','n_a_dwn','-------','n_b_up','-------','n_b_dwn'
      do i =1,d**2
       write(27,*) i,n_a_up(i),n_a_dwn(i),n_b_up(i),n_b_dwn(i)  ! orbitals: a and b
       write(28,*) i, real(phi_x(i)),real(phi_y(i)),real(phi_z(i)),real(zhi_x(i)),real(zhi_y(i)),real(zhi_z(i))
     
      enddo
!*********************************************************************************************   
     
     call get_mu
     f=0.0d0
	 do i=1,m_size
	 f=f+(1.0d0/(exp((w(i)-m_d)/temp)+1.0d0))
	 end do
     write(*,*)f,m_d

    count2 = 0
    sum2 = 0.0d0

    do i =1,m_size
    if(w(i).le.m_d) count2 =count2+1
    enddo
    write(60,*) count2
	
     
    !sum2 = 0.0d0
    !do i =1,d**2
    !sum2 = sum2 + n_a_up(i)+n_b_up(i)+n_a_dwn(i)+n_b_dwn(i)  !Total no.of particles!
    !enddo
    !write(*,*) sum2,'sum'
   

!! DOS calculation
 write(29,*) m_d
 !do i=1,4*d**2
    !write(29,*) i,w(i)
 !enddo 

double precision :: d_omega
d_omega=0.01d0

dos_en1 =  w(1)-m_d -0.5
do while (dos_en1.le. w(4*d**2)-m_d + 0.5)
dos = (0.0,0.0)
do i = 1, 4*d**2
dos = dos + (gama_d/pii)/(gama_d**2 + (dos_en1- (w(i) - m_d) )**2)
enddo
dos_en1 = dos_en1 + d_omega
write(100,*) dos_en1, real(dos) 
enddo

!!Projected DOS



dos_en =  w(1)-m_d -0.5
omega = 0
do while (dos_en.le. w(4*d**2)-m_d +0.5)
dos_aa_up = (0.0,0.0)
dos_aa_dwn = (0.0,0.0)
dos_bb_up = (0.0,0.0)
dos_bb_dwn = (0.0,0.0)
do j =1, d**2
do i = 1, 4*d**2

dos_aa_up = dos_aa_up+  H(2*j-1,i)*conjg(H(2*j-1,i))*((gama_d/pii)/(gama_d**2 + (dos_en- (w(i) - m_d)   )**2))
dos_aa_dwn = dos_aa_dwn+  H(2*d**2+2*j-1,i)*conjg(H(2*d**2+2*j-1,i))*((gama_d/pii)/(gama_d**2 + (dos_en- (w(i) - m_d)  )**2))


dos_bb_up = dos_bb_up + H(2*j,i)*conjg(H(2*j,i))*((gama_d/pii)/(gama_d**2 + (dos_en-  (w(i) - m_d) )**2))
dos_bb_dwn = dos_bb_dwn + H(2*d**2+2*j,i)*conjg(H(2*d**2+2*j,i))*((gama_d/pii)/(gama_d**2 + (dos_en-  (w(i) - m_d) )**2))



enddo
enddo
dos_en = dos_en + d_omega 
write(102,*) dos_en,  real(dos_aa_up),real(dos_bb_up)
write(103,*) dos_en, real(dos_aa_dwn),real(dos_bb_dwn)
write(104,*) dos_en, real(dos_aa_up+dos_aa_dwn+dos_bb_up+dos_bb_dwn)
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! IPR


 !l=0
 
 !BB2 = 0.0d0
 ipr2 = 0.0d0
 ipr_a = 0.0d0
 ipr_b = 0.0d0
 
  k =0
  
  en_mu_l = m_d -delta_mu
  en_mu_u = m_d +delta_mu

  do j =1,4*d**2

  if((w(j).le.en_mu_u).and.(w(j).ge.en_mu_l)) then
  k =k+1
  BB1 = 0.0d0
  BB2 = 0.0d0
  BB3 = 0.0d0
  BB4 = 0.0d0
  BB5 = 0.0d0
  BB6 = 0.0d0

  do i =1,2*d**2,2
  
  !write(*,*) 'yes'
   
 BB1  =BB1+((H(i,j) + H(i+1,j)+ H(i+2*d**2,j)+ H(i+(2*d**2)+1,j))*&
           conjg((H(i,j) + H(i+1,j)+ H(i+2*d**2,j)+ H(i+(2*d**2)+1,j))))**2
 BB2  =BB2+(H(i,j) + H(i+1,j)+ H(i+2*d**2,j)+ H(i+(2*d**2)+1,j))*&
           conjg((H(i,j) + H(i+1,j)+ H(i+2*d**2,j)+ H(i+(2*d**2)+1,j)))
 
 BB3 = BB3 + ((H(i,j)+ H(i+2*d**2,j))*conjg((H(i,j)+ H(i+2*d**2,j))))**2 
 BB4 = BB4 + (H(i,j)+ H(i+2*d**2,j))*conjg((H(i,j)+ H(i+2*d**2,j)))

 BB5 = BB5 + ((H(i+1,j)+ H(i+(2*d**2)+1,j))*conjg((H(i+1,j)+ H(i+(2*d**2)+1,j))))**2 
 BB6 = BB6 + (H(i+1,j)+ H(i+(2*d**2)+1,j))*conjg((H(i+1,j)+ H(i+(2*d**2)+1,j)))
 enddo
 ipr2=ipr2+(BB1/BB2**2)
 ipr_a = ipr_a + (BB3/BB4**2)
 ipr_b = ipr_b + (BB5/BB6**2)
 !write(*,*) ipr2,k
 endif

 !l=l+1
 !ipr(l) = BB1/BB2**2

 
 enddo
 !BB1 = 0.0d0
  !do i =1,d**2
     !BB1 = BB1+ipr(i)
     !write(31,*) i,ipr(i)/k
  !enddo
  write(31,*) V, ipr2/dble(k) !BB1/(BB2)**2
  write(41,*) V,ipr_a/dble(k) 
  write(42,*) V,ipr_b/dble(k) 
  !write(*,*) ipr3
  write(33,*) V,m_d,delta_mu,k
  write(33,*) '*****************************'
  do i=1,4*d**2
      write(33,*) i,w(i)
  enddo 




 ipr2 = 0.0d0
 ipr_a = 0.0d0
 ipr_b = 0.0d0
 
  k =0
  
  en_mu_l = m_d -delta_mu
  en_mu_u = m_d +delta_mu

  do j =1,4*d**2

  if((w(j).le.en_mu_u).and.(w(j).ge.en_mu_l)) then
  k =k+1
  BB1 = 0.0d0
  BB2 = 0.0d0
  BB3 = 0.0d0
  BB4 = 0.0d0
  BB5 = 0.0d0
  BB6 = 0.0d0

  do i =1,2*d**2,2
  
 
 BB1  =BB1+((H(i,j)*conjg(H(i,j))) + (H(i+1,j)*conjg(H(i+1,j))) &
          &+ ((H(i+2*d**2,j))*conjg(H(i+2*d**2,j))) +&
          & ((H(i+(2*d**2)+1,j))*conjg(H(i+(2*d**2)+1,j))))**2
 BB2  =BB2+(H(i,j)*conjg(H(i,j))) + (H(i+1,j)*conjg(H(i+1,j))) &
          &+ ((H(i+2*d**2,j))*conjg(H(i+2*d**2,j))) +&
          & ((H(i+(2*d**2)+1,j))*conjg(H(i+(2*d**2)+1,j)))

 BB3 = BB3 + ((H(i,j)*conjg(H(i,j))) + ((H(i+2*d**2,j))*conjg(H(i+2*d**2,j))))**2 
 BB4 = BB4 + (H(i,j)*conjg(H(i,j))) +&
           & ((H(i+(2*d**2),j))*conjg(H(i+(2*d**2),j)))

 BB5 = BB5 + ((H(i+1,j)*conjg(H(i+1,j)))+&
          & ((H(i+(2*d**2)+1,j))*conjg(H(i+(2*d**2)+1,j))))**2 

 BB6 = BB6 + (H(i+1,j)*conjg(H(i+1,j)))+&
          & ((H(i+(2*d**2)+1,j))*conjg(H(i+(2*d**2)+1,j)))
  
 enddo
 
 ipr2=ipr2+(BB1/BB2**2)
 ipr_a = ipr_a + (BB3/BB4**2)
 ipr_b = ipr_b + (BB5/BB6**2)

 endif

 enddo

  write(61,*) V, ipr2/dble(k) !BB1/(BB2)**2
   write(71,*) V,ipr_a/dble(k) 
  write(72,*) V,ipr_b/dble(k) 
   write(73,*) V,m_d,delta_mu,k
  

!!! Probability distribution 


 do i =1,d**2
      num(i) = n_a_up(i)+n_b_up(i)+n_a_dwn(i)+n_b_dwn(i)
      write(63,*) i,num(i)
 enddo
 
!gama = 0.050d0
prob(:) =0.0d0
numbr = 0.0d0
omega = 0
do while (numbr.le. 5.0d0)
!omega = omega+1
prob_sum = 0.0d0
do i = 1, d**2
prob_sum = prob_sum+ (gama/pii)/(gama**2 + (numbr-num(i))**2)
enddo
numbr = numbr +0.010d0 
write(101,*) numbr, prob_sum
enddo

!!!! Structure factors

 do i_nx = 1,d+1
      q_x = -pii+2.0d0*(i_nx-1)*pii/dble(d)
   do i_ny =1,d+1
      q_y =  -pii+2.0d0*(i_ny-1)*pii/dble(d)
   
  
      ssf = 0.0d0
      osf = 0.0d0
   do i =1,d**2
      do j =1,d**2
       arg = q_x*(pi(i)-pi(j))+q_y*(pj(i)-pj(j))
      !ssf = ssf + phi_z(i)*phi_z(j)*cmplx(cos(arg),sin(arg))
       ssf = ssf + (phi_x(i)*phi_x(j) +phi_y(i)*phi_y(j)+phi_z(i)*phi_z(j))*cmplx(cos(arg),sin(arg))
      !osf = osf + zhi_z(i)*zhi_z(j)*cmplx(cos(arg),sin(arg))
       osf = osf + (zhi_x(i)*zhi_x(j) +zhi_y(i)*zhi_y(j)+zhi_z(i)*zhi_z(j))*cmplx(cos(arg),sin(arg))
      enddo
   enddo
  write(85,*) q_x,q_y,real(ssf/(d**4))
  write(86,*) q_x,q_y,real(osf/(d**4))
   enddo
   enddo
    

!!!! Quantum Structure Factors

!Orbital

!do i_nx = 1,d+1
      !q_x = -pii+2.0d0*(i_nx-1)*pii/dble(d)
   !do i_ny =1,d+1
      !q_y =  -pii+2.0d0*(i_ny-1)*pii/dble(d)
   il = 0
 
      qosf = 0.0d0
     

 do i = 1,2*d**2,2
        il = il +1
        jl = 0
     do j = 1,2*d**2,2 
        jl = jl+1
      qof_xy1 = 0.0d0
      qof_xy2 = 0.0d0
      qof_xy3 = 0.0d0
      qof_xy4 = 0.0d0
      qof_xy5 = 0.0d0
      qof_xy6 = 0.0d0
      qof_xy7 = 0.0d0
      qof_xy8 = 0.0d0
      qof_z1 = 0.0d0
      qof_z2 = 0.0d0
      qof_z3 = 0.0d0
      qof_z4 = 0.0d0
      arg = q_x*(pi(il)-pi(jl))+q_y*(pj(il)-pj(jl))

        do nn = 1,4*d**2
            qof_xy1 = qof_xy1 + (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(i+1,nn))*H(i,nn)+ &
                    (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(i+(2*d**2)+1,nn))*H(i+(2*d**2),nn)
            qof_xy3 = qof_xy3 + (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(i+1,nn))*H(j+1,nn)+&
                    (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(i+((2*d**2)+1),nn))*H(j+1,nn)+&
                    (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(i+1,nn))*H(j+((2*d**2)+1),nn)+&
                    (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(i+((2*d**2)+1),nn))*H(j+((2*d**2)+1),nn)

            qof_xy5 = qof_xy5 + (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(i,nn))*H(i+1,nn)+ &
                    (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(i+(2*d**2)+1,nn))*H(i+2*d**2,nn)

            qof_xy7 = qof_xy7 + (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(i,nn))*H(j,nn)+&
                    (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(i+(2*d**2),nn))*H(j,nn)+&
                    (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(i,nn))*H(j+(2*d**2),nn)+&
                    (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(i+(2*d**2),nn))*H(j+(2*d**2),nn)

            qof_z1 = qof_z1 + (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(i,nn))*H(i,nn)+&
                     (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(i+(2*d**2),nn))*H(i+(2*d**2),nn)
            qof_z2 = qof_z2 + (1/(1+exp((w(n)-m_d)/temp)))*conjg(H(j,nn))*H(j,nn)+&
                     (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(j+(2*d**2),nn))*H(j+(2*d**2),nn)
            qof_z3 = qof_z3 + (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(i+1,nn))*H(i+1,nn)+&
                     (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(i+((2*d**2)+1),nn))*H(i+((2*d**2)+1),nn)
            
            qof_z4 = qof_z4 + (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(j+1,nn))*H(j+1,nn)+&
                     (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(j+((2*d**2)+1),nn))*H(j+((2*d**2)+1),nn)

             

            
        enddo
         do k = 1,4*d**2  
            qof_xy2 = qof_xy2 + (1/(1+exp((w(k)-m_d)/temp)))*conjg(H(j,k))*H(j+1,k)+&
                        (1/(1+exp((w(k)-m_d)/temp)))*conjg(H(j+(2*d**2),k))*H(j+(2*d**2)+1,k)
            qof_xy4 = qof_xy4 + (1-(1/(1+exp((w(k)-m_d)/temp))))*conjg(H(j,k))*H(i,k)+&
                    (1-(1/(1+exp((w(k)-m_d)/temp))))*conjg(H(j+(2*d**2),k))*H(i,k)+&
                    (1-(1/(1+exp((w(k)-m_d)/temp))))*conjg(H(j,k))*H(i+(2*d**2),k)+&
                    (1-(1/(1+exp((w(k)-m_d)/temp))))*conjg(H(j+(2*d**2),k))*H(i+(2*d**2),k)

            qof_xy6 = qof_xy6 + (1/(1+exp((w(k)-m_d)/temp)))*conjg(H(j+1,k))*H(j,k)+&
                        (1/(1+exp((w(k)-m_d)/temp)))*conjg(H(j+((2*d**2)+1),k))*H(j+(2*d**2),k)

            qof_xy8 = qof_xy8 + (1-(1/(1+exp((w(k)-m_d)/temp))))*conjg(H(j+1,k))*H(i+1,k)+&
                    (1-(1/(1+exp((w(k)-m_d)/temp))))*conjg(H(j+((2*d**2)+1),k))*H(i+1,k)+&
                    (1-(1/(1+exp((w(k)-m_d)/temp))))*conjg(H(j+1,k))*H(i+((2*d**2)+1),k)+&
                    (1-(1/(1+exp((w(k)-m_d)/temp))))*conjg(H(j+((2*d**2)+1),k))*H(i+((2*d**2)+1),k) 

           
         enddo
         
         !qosf = qosf + (0.50d0*(qof_xy1*qof_xy2 + qof_xy3*qof_xy4 + qof_xy5*qof_xy6 + qof_xy7*qof_xy8)+&
         !      0.250d0*(qof_z1*qof_z2 -qof_z3*qof_z2 -qof_z1*qof_z4 + qof_z3*qof_z4)) !*cmplx(cos(arg),sin(arg))
     
          qosf = (0.50d0*(qof_xy1*qof_xy2 + qof_xy3*qof_xy4 + qof_xy5*qof_xy6 + qof_xy7*qof_xy8)+&
               0.250d0*(qof_z1*qof_z2 -qof_z3*qof_z2 -qof_z1*qof_z4 + qof_z3*qof_z4)) 
     !
        !if (i.eq.j) then  
           !write(90,*) (0.50d0*(qof_xy1*qof_xy2 + qof_xy3*qof_xy4 + qof_xy5*qof_xy6 + qof_xy7*qof_xy8)+&
           !    0.250d0*(qof_z1*qof_z2 -qof_z3*qof_z2 -qof_z1*qof_z4 + qof_z3*qof_z4))
        !endif 
        !write (90,*) il,jl, real(qosf)  
 enddo
  
 enddo
         !write(185,*) q_x,q_y,real(qosf/(d**4)),imag(qosf/(d**4))
 !enddo
 !enddo
   

! Spin

      il = 0
 
      qsf = 0.0d0
     

 do i = 1,2*d**2,2
        il = il +1
        jl = 0
     do j = 1,2*d**2,2 
        jl = jl+1
        qsf_xy1 =0.0d0
        qsf_xy2 = 0.0d0
        qsf_xy3 = 0.0d0
        qsf_xy4 = 0.0d0
        qsf_z1 = 0.0d0
        qsf_z2 = 0.0d0
        qsf_z3 = 0.0d0
        qsf_z4 = 0.0d0
        qsf_z5 = 0.0d0
        qsf_z6 = 0.0d0
        qsf_z7 = 0.0d0
        qsf_z8 = 0.0d0
        qsf_z9 = 0.0d0
        qsf_z10 = 0.0d0
        qsf_z11 = 0.0d0
        qsf_z12 = 0.0d0
       
     do nn = 1,4*d**2
          
     qsf_xy1 = qsf_xy1 + (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(i+1,nn))*H(i+(2*d**2)+1,nn)+ &
                    (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(i+(2*d**2)+1,nn))*H(i+1,nn)
     qsf_xy2 = qsf_xy2 + (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(j,nn))*H(j+(2*d**2),nn)+ &
                    (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(j+(2*d**2),nn))*H(j,nn)
     
     qsf_xy3 = qsf_xy3 + (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(i,nn))*H(i+(2*d**2),nn)+ &
                    (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(i+(2*d**2),nn))*H(i,nn)

     qsf_xy4 = qsf_xy4 + (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(j+1,nn))*H(j+(2*d**2)+1,nn)+ &
                    (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(j+(2*d**2)+1,nn))*H(j+1,nn)

     qsf_z1 = qsf_z1 + (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(i,nn))*H(i,nn)+ &
                    (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(i+(2*d**2),nn))*H(i+(2*d**2),nn)

     qsf_z2 = qsf_z2 + (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(j,nn))*H(j,nn)+ &
                    (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(j+(2*d**2),nn))*H(j+(2*d**2),nn)
    
     qsf_z3 = qsf_z3 + (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(i+1,nn))*H(i+1,nn)
      
     qsf_z4 = qsf_z4 + (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(i+(2*d**2)+1,nn))*H(i+(2*d**2)+1,nn)

     qsf_z5 = qsf_z5 + (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(i,nn))*H(i,nn)
      
     qsf_z6 = qsf_z6 + (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(i+(2*d**2),nn))*H(i+(2*d**2),nn)

     qsf_z7 = qsf_z7 + (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(i+1,nn))*H(i+1,nn)+ &
                    (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(i+(2*d**2)+1,nn))*H(i+(2*d**2)+1,nn)

     qsf_z8 = qsf_z8 + (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(j,nn))*H(j,nn)+ &
                    (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(j+(2*d**2)+1,nn))*H(j+(2*d**2)+1,nn)

     qsf_z9 = qsf_z3 + (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(j+1,nn))*H(j+1,nn)
      
     qsf_z10 = qsf_z4 + (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(j+(2*d**2)+1,nn))*H(j+(2*d**2)+1,nn)

     qsf_z11 = qsf_z5 + (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(j,nn))*H(j,nn)
      
     qsf_z12 = qsf_z6 + (1/(1+exp((w(nn)-m_d)/temp)))*conjg(H(j+(2*d**2),nn))*H(j+(2*d**2),nn)

 

                    
    enddo
    qsf = 0.50d0*((qsf_xy1*qsf_xy2 + qsf_xy3*qsf_xy4)+(qsf_z1*qsf_z2 - qsf_z5*qsf_z10 -&
                        qsf_z6*qsf_z9 + qsf_z7*qsf_z8 -qsf_z3*qsf_z5 -qsf_z4*qsf_z12)) 

    ! write (990,*) il,jl, real(qsf)
  enddo
enddo

 call Create_SP_Density_Matrix
 call Calculate_two_point_correlations


do i_nx = 1,d+1
      q_x = -pii+2.0d0*(i_nx-1)*pii/dble(d)
   do i_ny =1,d+1
      q_y =  -pii+2.0d0*(i_ny-1)*pii/dble(d)
      Quant_orb_strfac = (0.0d0,0.0d0)
      Quant_spin_strfac =(0.0d0,0.0d0) 
 do i = 1,d**2
    do j =1,d**2
         arg = q_x*(pi(i)-pi(j))+q_y*(pj(i)-pj(j))
 Quant_orb_strfac =  Quant_orb_strfac + T_i_T_j(i,j)*cmplx(cos(arg),sin(arg))
 Quant_spin_strfac = Quant_spin_strfac +S_i_S_j(i,j)*cmplx(cos(arg),sin(arg))
    enddo
 enddo
  write(90,*) q_x,q_y,real(Quant_spin_strfac/(d**4)),imag(Quant_spin_strfac/(d**4))
             
  write(95,*) q_x,q_y, real(Quant_orb_strfac/(d**4)),imag(Quant_orb_strfac/(d**4))
    enddo
 enddo

   end   
!**********************************Main program ends******************************************


!*********************************************************************************************
!*****************Matrix construction + Diagonalization***************************************

    subroutine matrix_gen
    use mtmod
    use shared_arr
    use glob_var
	implicit none
    
	double precision m,free,newFE,esum1
	integer l,i,j,k,a,b,ii,ji,info,g,id,jd,flag
	integer i_i,j_j,s,ref
    double precision musum,mu,Fermi,get_mu,el_sum
   
  
  
    

    !Kinetic energy matrix

    t_arr(1)=t_aa
    t_arr(2)=t_ab
    t_arr(3)=t_ba
    t_arr(4)=t_bb


 
    k=1
	do i=1,shift
	 do j=1,shift
	  t(i,j)=t_arr(k)
	  k=k+1
	 end do
	end do


! Matrix Construction

     H=dcmplx(0.0d0,0.0d0)
    
	
     do l=1,d**2
      
         i=pi(l)
	     j=pj(l)	
	     ii=1
	     id=-1
	     ji=1
	     jd=-1
	     if (i.eq.1) id=-i+d
	     if (i.eq.d) ii=1-i
	     if (j.eq.1) jd=-j+d
	     if (j.eq.d) ji=1-j

!******* Off-diagonal blocks *************************

	   a=2*( (d*((i+ii)-1)+j )-1)+1
	   b=2*( (d*(i-1)+j )-1)+1

        do i_i=1,shift
	    do j_j=1,shift
	     H(a+(i_i-1),b+(j_j-1))=t(i_i,j_j)
        end do
	    end do
	   
	   a=2*( (d*((i+id)-1)+j )-1)+1
	   b=2*( (d*(i-1)+j )-1)+1

        do i_i=1,shift
	    do j_j=1,shift
	     H(a+(i_i-1),b+(j_j-1))=t(i_i,j_j)
        end do
	    end do
	   
	   a=2*((d*(i-1)+j)-1)+1
	   b=2*((d*(i-1)+j+ji)-1)+1	
       
        do i_i=1,shift
	    do j_j=1,shift
	     H(a+(i_i-1),b+(j_j-1))=t(i_i,j_j)
        end do
	    end do

	   a=2*((d*(i-1)+j)-1)+1
	   b=2*((d*(i-1)+j+jd)-1)+1

        do i_i=1,shift
	    do j_j=1,shift
	     H(a+(i_i-1),b+(j_j-1))=t(i_i,j_j)
         end do
	     end do

 !************ Diagonal block************************
      
         a=2*((d*(i-1)+j)-1)+1
	     b=2*((d*(i-1)+j)-1)+1
          
         
       s=0
       call mfmat(l)
       call rannum(r)
	   do i_i=1,shift
	   do j_j=1,shift
        s=s+1
        if (i_i.eq.j_j) then 
       
         !e_onsite = r*V 
         !write(*,*) e_onsite   
        H(a+(i_i-1),b+(j_j-1))=up_up(l,s)+dis(l)
        H(a+(i_i-1)+2*d**2,b+(j_j-1)+2*d**2) = down_down(l,s)+dis(l)
        else 
        H(a+(i_i-1),b+(j_j-1))=up_up(l,s)
        H(a+(i_i-1)+2*d**2,b+(j_j-1)+2*d**2) =down_down(l,s)
        endif
        H(a+(i_i-1),b+(j_j-1)+2*d**2)=up_down(l,s)
        H(a+(i_i-1)+2*d**2,b+(j_j-1)) =down_up(l,s)
       enddo
       enddo	
       
        
    enddo

       !do i =1,m_size
          !do j=1,m_size
             !write(30,*) i,j,real(H(i,j))
          !enddo
       !enddo
       
!Diagonalization
      
        call zheev ('V','U',4*d**2,H,4*d**2,w,work,2*(4*d**2)-1,rwork,info)
        if (info.ne.0) print*,"info=",info
		!print*,real(w)


 return
 end

!*********************************************************************************************

   
!********************************Chemical Potential*******************************************

    subroutine get_mu
	use shared_arr
	use glob_var
	implicit none
	double precision f0, f, fL2, fR, m, mR, mL, rtmp
	integer i,omega
    !write(*,*) 'passed1',  U,fill
    
     omega = 0
	 mR = maxval(w)       !right-side chemical potential
    
	fR=0.0d0
	do i=1,m_size
   
	fR=fR+(1.0d0/(exp((w(i)-mR)/temp)+1.0d0))
	end do
    
	 mL = minval(w)       !left-side chemical potential
     
	 fL2=0.0d0
 	do i=1,m_size
 	fL2=fL2+(1.0d0/(exp((w(i)-mL)/temp)+1.0d0))
	end do
    
	m_d = 0.5d0*(mL+mR)    !middle chemical potential
    
	f=0.0d0
	do i=1,m_size
	f=f+(1.0d0/(exp((w(i)-m_d)/temp)+1.0d0))
	end do
	
	do while(abs((f-fill)/dble(d**2)).ge.1e-4)
    omega=omega+1
     !write(*,*) (f-fill), '(f-fill)', omega,'omega'
	 m_d = 0.5d0*(mL+mR)
    
	f=0.0d0
	do i=1,m_size
	f=f+(1.0d0/(exp((w(i)-m_d)/temp)+1.0d0))
	end do	  
	 if(f.gt.fill)then
	  !if middle filling is above target, make it the new right bound.
	  mR = m_d
	  fR = f
	 elseif(f.lt.fill)then
	  !if middle filling is below target, make it the new left bound. 
	  mL = m_d
	  fR = f
	 endif
	enddo
   
    !write (*,*) fill, m_d, "mu value"
	
	return
	end

!*********************************************************************************************


!*********************MF parameters in 2*2 matrices*******************************************

    subroutine mfmat(ref)
    use shared_arr
    use glob_var

 
 integer i,l,j,ref
 up_up(ref,1) = -(A1*phi_z(ref)+B1*zhi_z(ref))
 !up_up(ref,2) = -(B1/2)*cmplx(zhi_x(ref),-zhi_y(ref))
 up_up(ref,2) = -(B1)*(zhi_x(ref)-i_complex*zhi_y(ref))
 !up_up(ref,3) = -(B1/2)*cmplx(zhi_x(ref),zhi_y(ref))
 up_up(ref,3) = -(B1)*(zhi_x(ref) +i_complex*zhi_y(ref))
 up_up(ref,4) =  (-A1*phi_z(ref)+B1*zhi_z(ref))

 !up_down(ref,1) = -(A1/2)*cmplx(phi_x(ref),-phi_y(ref))
 up_down(ref,1) = -(A1)*(phi_x(ref)-i_complex*phi_y(ref))
 up_down(ref,2) = 0.0d0
 up_down(ref,3) = 0.0d0
 !up_down(ref,4) =-(A1/2)*cmplx(phi_x(ref),-phi_y(ref))
  up_down(ref,4) =-(A1)*(phi_x(ref)-i_complex*phi_y(ref))

 !down_up(ref,1) =-(A1/2)*cmplx(phi_x(ref),phi_y(ref))
  down_up(ref,1) =-(A1)*(phi_x(ref) + i_complex*phi_y(ref))
 down_up(ref,2) = 0.0d0
 down_up(ref,3) = 0.0d0
 !down_up(ref,4) =-(A1/2)*cmplx(phi_x(ref),phi_y(ref))
  down_up(ref,4) =-(A1)*(phi_x(ref)+i_complex*phi_y(ref))

 down_down(ref,1)= (A1*phi_z(ref)-B1*zhi_z(ref))
 !down_down(ref,2)= -(B1/2)*cmplx(zhi_x(ref),-zhi_y(ref))
  down_down(ref,2)= -(B1)*(zhi_x(ref)-i_complex*zhi_y(ref))
 !down_down(ref,3)= -(B1/2)*cmplx(zhi_x(ref),zhi_y(ref))
  down_down(ref,3)= -(B1)*(zhi_x(ref)+i_complex*zhi_y(ref))
 down_down(ref,4)= (A1*phi_z(ref)+B1*zhi_z(ref))

 return
 end

!*********************************************************************************************


!****************************Calculation of MF parameters*************************************

 subroutine mf
 use shared_arr
 use glob_var
 integer ::i,j,kk,l
 complex*16::phi_sum1,zhi_sum1, phi_sum2,zhi_sum2
 complex*16::phi_sum3,zhi_sum3
 
 
  !write(*,*) 'passed2', m_d,'mu'
 
 AA1 = 0.0d0
 AA2 = 0.0d0
 AA3 = 0.0d0
 AA4 = 0.0d0
 AA5 = 0.0d0
 AA6 = 0.0d0
 AA7 = 0.0d0
 AA8 = 0.0d0
 AA9 = 0.0d0
 AA10 = 0.0d0
 AA11= 0.0d0
 AA12 = 0.0d0
 
 
 l=0
 
 do i =1,2*d**2,2
    !write(*,*) i,i+2*d**2,i+(2*d**2)+1
  AA1 = 0.0d0
  AA2 = 0.0d0
  AA3 = 0.0d0
  AA4 = 0.0d0
  AA5 = 0.0d0
  AA6 = 0.0d0
  AA7 = 0.0d0
  AA8 = 0.0d0
  AA9 = 0.0d0
  AA10 = 0.0d0
  AA11= 0.0d0
  AA12 = 0.0d0
  phi_sum1 = 0.0d0
  phi_sum2 = 0.0d0
  phi_sum3 = 0.0d0
  zhi_sum1 = 0.0d0
  zhi_sum2 = 0.0d0
  zhi_sum3 = 0.0d0



 do j =1,m_size
 AA1  = AA1+H(i,j)*(1/(1+exp((w(j)-m_d)/temp)))*conjg (H(i,j)) !n_a_up
 AA2  = AA2+H(i+2*d**2,j)*(1/(1+exp((w(j)-m_d)/temp)))*conjg (H(i+2*d**2,j)) !n_a_dwn
 AA3  = AA3+H(i+1,j)*(1/(1+exp((w(j)-m_d)/temp)))*conjg (H(i+1,j)) !n_b_up
 AA4  = AA4+H(i+(2*d**2)+1,j)*(1/(1+exp((w(j)-m_d)/temp)))*conjg (H(i+(2*d**2)+1,j))!n_b_dwn
 AA5  = AA5+H(i,j)*(1/(1+exp((w(j)-m_d)/temp)))*conjg(H(i+2*d**2,j))
 AA6  = AA6+H(i+2*d**2,j)*(1/(1+exp((w(j)-m_d)/temp)))*conjg (H(i,j))
 AA7  = AA7+H(i+1,j)*(1/(1+exp((w(j)-m_d)/temp)))*conjg (H(i+(2*d**2)+1,j))
 AA8  = AA8+H(i+(2*d**2)+1,j)*(1/(1+exp((w(j)-m_d)/temp)))*conjg (H(i+1,j)) 
 AA9  = AA9+H(i,j)*(1/(1+exp((w(j)-m_d)/temp)))*conjg (H(i+1,j)) 
 AA10 = AA10+H(i+1,j)*(1/(1+exp((w(j)-m_d)/temp)))*conjg (H(i,j)) 
 AA11 = AA11+H(i+2*d**2,j)*(1/(1+exp((w(j)-m_d)/temp)))*conjg (H(i+(2*d**2)+1,j)) 
 AA12 = AA12+H(i+(2*d**2)+1,j)*(1/(1+exp((w(j)-m_d)/temp)))*conjg (H(i+2*d**2,j)) 
 enddo


 
 phi_sum1 = 0.50d0*(AA5+AA6+AA7+AA8)
 phi_sum2 = 0.50d0*i_complex*(-AA5+AA6-AA7+AA8)  
 phi_sum3 = 0.50d0*(AA1+AA3-AA2-AA4)  !actual MF parameters
 zhi_sum1 = 0.50d0*(AA9+AA10+AA11+AA12)
 zhi_sum2 = 0.50d0*i_complex*(-AA9+AA10-AA11+AA12)
 zhi_sum3 = 0.50d0*(AA1+AA2-AA3-AA4) 
 l=l+1
  write(171,*) l,AA5,AA6,AA7,AA8
   n_a_up(l) = AA1
   n_a_dwn(l) = AA2
   n_b_up(l) = AA3
   n_b_dwn(l) = AA4
   
   !write(*,*) l,AA1,AA2 !
   !write(*,*) l,AA3,AA4
 
 phi_x(l) = phi_sum1
 phi_y(l) = phi_sum2
 phi_z(l) = phi_sum3
 zhi_x(l) = zhi_sum1
 zhi_y(l) = zhi_sum2
 zhi_z(l) = zhi_sum3
 !write(*,*) phi_sum2, phi_sum1,zhi_sum2,zhi_sum3

 enddo
 !pause
 return
 end

!*********************************************************************************************

subroutine Create_SP_Density_Matrix
use shared_arr
use glob_var
implicit none


integer::alpha1,beta1,n1

do alpha1=1,4*d**2
do beta1=1,4*d**2

SP_Density_Matrix(alpha1,beta1)=(0.0d0,0.0d0)
do n1=1,4*d**2
SP_Density_Matrix(alpha1,beta1) = SP_Density_Matrix(alpha1,beta1) + conjg(H(alpha1,n1))*H(beta1,n1)* &
					(1.0d0/ ( exp(   (w(n1)-m_d)*(1.0d0/temp)  ) + 1.0d0  ) )
enddo

enddo
enddo

return
end



!******************************************************************************************

complex*16 function Two_Particle_Density_Matrix(alpha_pass, beta_pass, gamma_pass, delta_pass)
!use mtmod
use shared_arr
use glob_var
implicit none
complex*16 :: delta_gamma_beta
integer :: alpha_pass, beta_pass, gamma_pass, delta_pass

if (gamma_pass .eq. beta_pass) then
delta_gamma_beta = (1.0d0,0.0d0)
else
delta_gamma_beta = (0.0d0,0.0d0)
endif

Two_Particle_Density_Matrix =  ( SP_Density_Matrix(alpha_pass,beta_pass)*SP_Density_Matrix(gamma_pass,delta_pass) ) + &
		 ( (SP_Density_Matrix(alpha_pass,delta_pass)) * ( delta_gamma_beta - SP_Density_Matrix(gamma_pass,beta_pass ) ) )


return 
end function Two_Particle_Density_Matrix

!******************************************************************************************

!*******************************************************************************************

subroutine Calculate_two_point_correlations
use shared_arr
use glob_var
implicit none

integer :: i,j,i_row, i_col, j_row, j_col
integer :: i_row_orb, i_col_orb, j_row_orb, j_col_orb
integer :: i_row_spin, i_col_spin, j_row_spin, j_col_spin
integer :: i1,i2,j1,j2

complex*16 :: Two_Particle_Density_Matrix
complex*16 :: temp_val_Sz_Sz, temp_val_Sx_Sx, temp_val_Sy_Sy
complex*16 :: temp_val_Tz_Tz, temp_val_Tx_Tx, temp_val_Ty_Ty
!BASIS USED FOR SINGLE SITE OPR's
!{a,up  b,up  a,dn  b,dn}
!{1  2  3   4}
!orbital a = 1, b=2
!spin up =1, dn =2
!onsite_basis_index = orb + (spin-1)*2

!basis_index = orb + ( site -1 )*2 + (spin -1)*2*d**2
!site = iy + (ix-1)*d [connected with calculation of pi,pj above!!]
  Sz_local(:,:) = 0.0d0
  Sx_local(:,:) = 0.0d0
  Sy_local(:,:) = 0.0d0

  Tz_local(:,:) = 0.0d0
  Tx_local(:,:) = 0.0d0
  Ty_local(:,:) = 0.0d0

 Sz_local(1,1) = (0.50d0,0.0d0)
 Sz_local(2,2) = (0.50d0,0.0d0)
 Sz_local(3,3) = (-0.50d0,0.0d0)
 Sz_local(4,4) = (-0.50d0,0.0d0)

 Sx_local(1,3) = (0.50d0,0.0d0) 
 Sx_local(2,4) = (0.50d0,0.0d0) 
 Sx_local(3,1) = (0.50d0,0.0d0) 
 Sx_local(4,2) = (0.50d0,0.0d0) 

 Sy_local(1,3) = (0.0d0,-0.50d0) 
 Sy_local(2,4) = (0.0d0,-0.50d0) 
 Sy_local(3,1) = (0.0d0,0.50d0) 
 Sy_local(4,2) = (0.0d0,0.50d0)

 Tz_local(1,1) = (0.50d0,0.0d0)
 Tz_local(2,2) = (-0.50d0,0.0d0)
 Tz_local(3,3) = (0.50d0,0.0d0)
 Tz_local(4,4) = (-0.50d0,0.0d0)

 Tx_local(1,2) = (0.50d0,0.0d0) 
 Tx_local(2,1) = (0.50d0,0.0d0) 
 Tx_local(3,4) = (0.50d0,0.0d0) 
 Tx_local(4,3) = (0.50d0,0.0d0) 

 Ty_local(1,2) = (0.0d0,-0.50d0) 
 Ty_local(2,1) = (0.0d0,0.50d0) 
 Ty_local(3,4) = (0.0d0,-0.50d0) 
 Ty_local(4,3) = (0.0d0,0.50d0)


 S_i_S_j(:,:) = (0.0d0,0.0d0) 
 T_i_T_j(:,:) = (0.0d0,0.0d0) 

 do i =1,d**2   
     do j =1,d**2
     temp_val_Sz_Sz = (0.0d0,0.0d0) 
	 temp_val_Sx_Sx = (0.0d0,0.0d0) 
 	 temp_val_Sy_Sy = (0.0d0,0.0d0) 
	 temp_val_Tz_Tz = (0.0d0,0.0d0) 
	 temp_val_Tx_Tx = (0.0d0,0.0d0) 
 	 temp_val_Ty_Ty = (0.0d0,0.0d0) 

	
	
	do i_row_orb=1,2
		do i_row_spin=1,2
			i_row = i_row_orb + (i_row_spin - 1)*2

		do i_col_orb=1,2
		do i_col_spin=1,2
			i_col = i_col_orb + (i_col_spin - 1)*2


	do j_row_orb=1,2
		do j_row_spin=1,2
			j_row = j_row_orb + (j_row_spin - 1)*2

		do j_col_orb=1,2
		do j_col_spin=1,2
			j_col = j_col_orb + (j_col_spin - 1)*2


		i1=i_row_orb + ( i -1 )*2 + (i_row_spin -1)*2*d**2
		i2=i_col_orb + ( i -1 )*2 + (i_col_spin -1)*2*d**2
		j1=j_row_orb + ( j -1 )*2 + (j_row_spin -1)*2*d**2
		j2=j_col_orb + ( j -1 )*2 + (j_col_spin -1)*2*d**2


		!write (*,*) Two_Particle_Density_Matrix(i1,i2,j1,j2)

		if ( (Sz_local(i_row,i_col) .ne. (0.0d0,0.0d0))  .and.  (Sz_local(j_row,j_col) .ne. (0.0d0,0.0d0)) ) then
		temp_val_Sz_Sz = temp_val_Sz_Sz + Sz_local(i_row,i_col)*Sz_local(j_row,j_col)*Two_Particle_Density_Matrix(i1,i2,j1,j2)	
		endif

        if ( (Sx_local(i_row,i_col) .ne. (0.0d0,0.0d0))  .and.  (Sx_local(j_row,j_col) .ne. (0.0d0,0.0d0)) ) then
		temp_val_Sx_Sx = temp_val_Sx_Sx + Sx_local(i_row,i_col)*Sx_local(j_row,j_col)*Two_Particle_Density_Matrix(i1,i2,j1,j2)		
		endif

        if ( (Sy_local(i_row,i_col) .ne. (0.0d0,0.0d0))  .and.  (Sy_local(j_row,j_col) .ne. (0.0d0,0.0d0)) ) then
		temp_val_Sy_Sy = temp_val_Sy_Sy + Sy_local(i_row,i_col)*Sy_local(j_row,j_col)*Two_Particle_Density_Matrix(i1,i2,j1,j2)		
		endif

		if ( (Tz_local(i_row,i_col) .ne. (0.0d0,0.0d0))  .and.  (Tz_local(j_row,j_col) .ne. (0.0d0,0.0d0)) ) then
		temp_val_Tz_Tz = temp_val_Tz_Tz + Tz_local(i_row,i_col)*Tz_local(j_row,j_col)*Two_Particle_Density_Matrix(i1,i2,j1,j2)	
		endif

        if ( (Tx_local(i_row,i_col) .ne. (0.0d0,0.0d0))  .and.  (Tx_local(j_row,j_col) .ne. (0.0d0,0.0d0)) ) then
		temp_val_Tx_Tx = temp_val_tx_Tx + Tx_local(i_row,i_col)*Tx_local(j_row,j_col)*Two_Particle_Density_Matrix(i1,i2,j1,j2)		
		endif

        if ( (Ty_local(i_row,i_col) .ne. (0.0d0,0.0d0))  .and.  (Ty_local(j_row,j_col) .ne. (0.0d0,0.0d0)) ) then
		temp_val_Ty_Ty = temp_val_Ty_Ty + Ty_local(i_row,i_col)*Ty_local(j_row,j_col)*Two_Particle_Density_Matrix(i1,i2,j1,j2)		
		endif
		


	enddo
	enddo


	enddo
	enddo	
   
 enddo
 enddo
 
enddo
enddo

 write(200,*) i,j,(temp_val_Sy_Sy+temp_val_Sx_Sx+temp_val_Sz_Sz) , (temp_val_Ty_Ty+temp_val_Tx_Tx+temp_val_Tz_Tz) 
 
  S_i_S_j(i,j) = (temp_val_Sy_Sy+temp_val_Sx_Sx+temp_val_Sz_Sz) 
  T_i_T_j(i,j) = (temp_val_Ty_Ty+temp_val_Tx_Tx+temp_val_Tz_Tz)  
enddo
enddo

return
end

!*******************************************************************************************


!*********************************************************************************************
!Random no. generator
!*********************************************************************************************
    	subroutine rannum(r1)
    	use mtmod
    	implicit none
    	 double precision r1
    	 r1=grnd()
     	return
     	end
