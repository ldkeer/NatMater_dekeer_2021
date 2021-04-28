      ! Configurations + identify types of atoms
      
      !Create ConnectivityMatrix
      allocate (ConnectivityMatrix(counter_ConnectivityArray,8))
      ConnectivityMatrix=0

      do i=1,counter_ConnectivityArray,1
              
          counter=0
          counter2=0
          do while (counter.lt.4)
              
              do j=1,4*counter_ConnectivityArray,1
                  
                  if (ConnectivityArray_left(j).eq.i) then
                      
                      
                      if (ConnectivityArray_left(j).eq.ConnectivityArray_right(j)) then
                          counter2=1
                      endif
                      
                      counter=counter+1
                      ConnectivityMatrix(i,counter)=ConnectivityArray_right(j)
                      ConnectivityMatrix(i,counter+4)=j
                      
                  endif
                  
                  if (ConnectivityArray_right(j).eq.i) then
                      
                      if (counter2.eq.1) then
                          counter=counter+1
                          ConnectivityMatrix(i,counter)=0
                          ConnectivityMatrix(i,counter+4)=0
                          counter2=0
                      else
                      
                      counter=counter+1
                      ConnectivityMatrix(i,counter)=ConnectivityArray_left(j)
                      ConnectivityMatrix(i,counter+4)=j
                      
                      endif
                      
                  endif
                  
              enddo
              
          enddo
          
      enddo
                  
      !Depth-first strategy to identify molecules
      !we go through connectivity matrix molecule per molecule and create copy matrix/arrays to identify all the atom(s) (group(s))
      !after identifying, we create internal coordinates per molecule (+adaption of angles in case cycles are formed afterwards)
      !path array defines the current sequence of CPs in a molecule
      !every time we add a new CPs to the path array, we first turn it's loose segments (not connected to another CP) into their internal coordinates

      allocate (Unimolecular_molecule(imax))  !to make sure we have all intramolecular paths for a certain molecule only once
      allocate (unimolecular_position(imax))  !to store number of rows in unimolecular matrix for every intramolecular path
      allocate (path(imax))   !defines current path in a molecule
      molecules_counter=1
      path=0
      counter=0
      unimolecular_counter=0
      unimolecular_position=0
      unimolecular_molecule=0
      
      allocate (Configuration(imax,5))    !column 1->type;column 2->number of atom considered;column 3/4/5->connectivities
      Configuration=0
      Configuration(1,1)=4 !type
      Configuration(1,2)=1
      counter_rowConfiguration=2
      
      allocate (Topology_copy(imax,jmax)) !to translate information to 'atom (group) numbers'
      allocate (ConnectivityArray_left_copy(imax))
      allocate (ConnectivityArray_right_copy(imax))
      Topology_copy=0
      ConnectivityArray_left_copy=0
      ConnectivityArray_right_copy=0
      
      do i=1,4*counter_ConnectivityArray,1
          if (ConnectivityArray_left(i).eq.1) then
              ConnectivityArray_left_copy(i)=1
          endif
          if (ConnectivityArray_right(i).eq.1) then
              ConnectivityArray_right_copy(i)=1
          endif
      enddo
      
      counter_visualization=1
      
      i=1
      j=1
      
      do while (molecules_counter.le.amountofMolecules)   !until all molecules are identified
          
          if (((ConnectivityMatrix(i,1)+ConnectivityMatrix(i,2)+ConnectivityMatrix(i,3)+ConnectivityMatrix(i,4)).eq.0)
     &.AND.(ConnectivityMatrix(i,5).ne.0).AND.(ConnectivityMatrix(i,6).ne.0).AND.(ConnectivityMatrix(i,7).ne.0)
     &.AND.(ConnectivityMatrix(i,8).ne.0)) then  !case 1: 4 'loose' segments; 1 molecule
!             
              do j=1,4,1 !'4 segments'
                  
                  if (counter_rowConfiguration.eq.2) then
                      Configuration(2,1)=Topology(ConnectivityMatrix(i,j+4),1)
                      Configuration(2,2)=2
                      Configuration(2,3)=1
                  else
                  if (counter_rowConfiguration.eq.3) then
                      Configuration(3,1)=Topology(ConnectivityMatrix(i,j+4),1)
                      Configuration(3,2)=3
                      Configuration(3,3)=1
                      Configuration(3,4)=2
                  else
                      Configuration(counter_rowConfiguration,1)=Topology(ConnectivityMatrix(i,j+4),1)
                      Configuration(counter_rowConfiguration,2)=4
                      Configuration(counter_rowConfiguration,3)=1
                      Configuration(counter_rowConfiguration,4)=2
                      Configuration(counter_rowConfiguration,5)=3

                  endif
                  endif
                      counter_rowConfiguration=counter_rowConfiguration+1
              enddo
              
              do k=1,4,1
                  ConnectivityMatrix(i,k+4)=0
              enddo
              
              if (molecules_counter.eq.(1+(counter_visualization-1)*ceiling(dble(amountofmolecules)/dble(molecules_output)))) then
              
!              write (5000+molecules_counter,59) Transpose(Configuration)
              allocate (Configuration_output(counter_rowConfiguration-1,9))   !input for gaussview; c1->bond length;c2->0/1;c3->bond angle;c4->0/1;c5->dihedral angle;c6->0/1;c7/8/9->connectivities 
              Configuration_output=0
              do i=1,counter_rowConfiguration-1,1
                  Configuration_output(i,7)=Configuration(i,3)
                  Configuration_output(i,8)=Configuration(i,4)
                  Configuration_output(i,9)=Configuration(i,5)
                  if (i.gt.1) then
                      include 'SiO.for'
                      Configuration_output(i,1)=muselected2
                      Configuration_output(i,2)=1
                  endif
                  if (i.gt.2) then
                  if ((Configuration(Configuration(i,2),1).eq.4).AND.(Configuration(Configuration(i,3),1).eq.3)
     &.AND.(Configuration(Configuration(i,4),1).eq.4)) then
                      include 'SiOSi.for'
                      Configuration_output(i,3)=muselected2
                  else
                      include 'OSiO.for'
                      Configuration_output(i,3)=muselected2
                  endif
                      Configuration_output(i,4)=1
                  endif                  
                  if (i.gt.3) then
                      include 'SiOSiO.for'
                      Configuration_output(i,5)=muselected2
                      Configuration_output(i,6)=1
                  endif     

              enddo
              write (6000+100*PopulationTimer3+counter_visualization,60) Transpose(Configuration_output)
              deallocate (Configuration_output)
              allocate (Configuration_type(counter_rowConfiguration-1))   !atom(group) type (to combine with Configuration_output as input Gaussview)
              do i=1,counter_rowConfiguration-1,1
                  if (Configuration(i,1).eq.1) then
                      Configuration_type(i)='B'
                  else
                  if (Configuration(i,1).eq.2) then
                      Configuration_type(i)='C'
                  else
                  if (Configuration(i,1).eq.3) then
                      Configuration_type(i)='O'
                  else
                      Configuration_type(i)='Si'
                  endif
                  endif
                  endif
              enddo
              write (7000+100*PopulationTimer3+counter_visualization,61) Configuration_type
              deallocate (Configuration_type)
              
              !can never have an intramolecular reaction so no adaption of configuration_output needed
              
              counter_visualization=counter_visualization+1
              endif
                  
              molecules_counter=molecules_counter+1
              Configuration=0
              Configuration(1,1)=4
              Configuration(1,2)=1
              counter_rowConfiguration=2
              ConnectivityArray_left_copy=0
              ConnectivityArray_right_copy=0
              Topology_copy=0
              if (molecules_counter.gt.amountofmolecules) then
              else
              !Add 'new first Si atom' to copy matrix/arrays)
              StopCrit3=0
              i=0
              do while (StopCrit3.eq.0)
                  i=i+1
                  if ((ConnectivityMatrix(i,5)+ConnectivityMatrix(i,6)+ConnectivityMatrix(i,7)+ConnectivityMatrix(i,8)).ne.0) then
                      StopCrit3=1
                  endif
              enddo
              do j=1,4*counter_ConnectivityArray,1
                  if (ConnectivityArray_left(j).eq.i) then
                      ConnectivityArray_left_copy(j)=1
                  endif
                  if (ConnectivityArray_right(j).eq.i) then
                      ConnectivityArray_right_copy(j)=1
                  endif
              enddo              
              j=1
              endif
              
          else
          if ((ConnectivityMatrix(i,j).eq.0).AND.(ConnectivityMatrix(i,j+4).ne.0)) then !case 2: 'free segment', always attached to Si that is already added

                  if (counter_rowConfiguration.eq.2) then
                      Configuration(2,1)=Topology(ConnectivityMatrix(i,j+4),1)
                      Configuration(2,2)=2
                      Configuration(2,3)=1
                      Topology_copy(ConnectivityMatrix(i,j+4),1)=2
                  else
                  if (counter_rowConfiguration.eq.3) then
                      Configuration(3,1)=Topology(ConnectivityMatrix(i,j+4),1)
                      Configuration(3,2)=3
                      Configuration(3,3)=1
                      Configuration(3,4)=2
                      Topology_copy(ConnectivityMatrix(i,j+4),1)=3
                  else
                      Configuration(counter_rowConfiguration,1)=Topology(ConnectivityMatrix(i,j+4),1)
                      Topology_copy(ConnectivityMatrix(i,j+4),1)=counter_rowConfiguration
                      Configuration(counter_rowConfiguration,2)=counter_rowConfiguration
                      Configuration(counter_rowConfiguration,3)=ConnectivityArray_left_copy(ConnectivityMatrix(i,j+4))
                      
                      if (ConnectivityArray_left_copy(ConnectivityMatrix(i,j+4)).eq.1) then
                          Configuration(counter_rowConfiguration,4)=2
                          Configuration(counter_rowConfiguration,5)=3
                      else
                      if (ConnectivityArray_left_copy(ConnectivityMatrix(i,j+4)).eq.2) then
                          Configuration(counter_rowConfiguration,4)=1
                          Configuration(counter_rowConfiguration,5)=3
                      else
                          Configuration(counter_rowConfiguration,4)
     &=Configuration(ConnectivityArray_left_copy(ConnectivityMatrix(i,j+4)),3)
                          Configuration(counter_rowConfiguration,5)
     &=Configuration(ConnectivityArray_left_copy(ConnectivityMatrix(i,j+4)),4)
                      endif
                      endif
                      
                  endif
                  endif

                  counter_rowConfiguration=counter_rowConfiguration+1
                  ConnectivityMatrix(i,j+4)=0
                  
                  if (j.lt.4) then
                      j=j+1   !1 to the right in connectivitymatrix
                  else
                      if (((ConnectivityMatrix(i,1)+ConnectivityMatrix(i,2)+ConnectivityMatrix(i,3)+ConnectivityMatrix(i,4)).eq.0)
     &.AND.((ConnectivityMatrix(i,5)+ConnectivityMatrix(i,6)+ConnectivityMatrix(i,7)+ConnectivityMatrix(i,8)).eq.0)) then !no further connectivities
                          if (counter.eq.0) then  !this should never happen I think, because from the moment you have a CP, counter becomes 1 and a 'free segment' can never be on its own 
                                  if (molecules_counter.eq.(1+(counter_visualization-1)*ceiling(dble(amountofmolecules)
     &/dble(molecules_output)))) then
                           !       write (5000+molecules_counter,59) Transpose(Configuration)
                                  allocate (Configuration_output(counter_rowConfiguration-1,9))
                                  Configuration_output=0
                                  do i=1,counter_rowConfiguration-1,1
                                      Configuration_output(i,7)=Configuration(i,3)
                                      Configuration_output(i,8)=Configuration(i,4)
                                      Configuration_output(i,9)=Configuration(i,5)
                                  if (i.gt.1) then
                                      include 'SiO.for'
                                      Configuration_output(i,1)=muselected2
                                      Configuration_output(i,2)=1
                                  endif
                                  if (i.gt.2) then
                                  if ((Configuration(Configuration(i,2),1).eq.4).AND.(Configuration(Configuration(i,3),1).eq.3)
     &.AND.(Configuration(Configuration(i,4),1).eq.4)) then
                                      include 'SiOSi.for'
                                      Configuration_output(i,3)=muselected2
                                  else
                                      include 'OSiO.for'
                                      Configuration_output(i,3)=muselected2
                                  endif
                                      Configuration_output(i,4)=1
                                  endif                  
                                  if (i.gt.3) then
                                      include 'SiOSiO.for'
                                      Configuration_output(i,5)=muselected2
                                      Configuration_output(i,6)=1
                                  endif     

                                  enddo
                                  
                                  !Adaption of configuration in case of intramolecular reactions
                                  if (unimolecular_counter.ne.0) then
                                  !    include 'unimolecular.for'
                                  endif
                                  unimolecular_counter=0  !amount of intramolecular paths (loops) per molecule
                                  unimolecular_position=0 !nr of rows in Unimolecular matrix representing intramolecular paths per molecule
                                  unimolecular_molecule=0 !to make sure we have all intramolecular paths for a certain molecule only once 
                                  
                                  write (6000+100*PopulationTimer3+counter_visualization,60) Transpose(Configuration_output)
                                  deallocate (Configuration_output)
                                  allocate (Configuration_type(counter_rowConfiguration-1))
                                  do i=1,counter_rowConfiguration-1,1
                                      if (Configuration(i,1).eq.1) then
                                          Configuration_type(i)='B'
                                      else
                                      if (Configuration(i,1).eq.2) then
                                          Configuration_type(i)='C'
                                      else
                                      if (Configuration(i,1).eq.3) then
                                          Configuration_type(i)='O'
                                      else
                                          Configuration_type(i)='Si'
                                      endif
                                      endif
                                      endif
                                  enddo
                                  write (7000+100*PopulationTimer3+counter_visualization,61) Configuration_type
                                  deallocate (Configuration_type)
                                  
                                  counter_visualization=counter_visualization+1
                                  endif
                                  
                                  molecules_counter=molecules_counter+1
                                  Configuration=0
                                  Configuration(1,1)=4
                                  Configuration(1,2)=1
                                  counter_rowConfiguration=2
              ConnectivityArray_left_copy=0
              ConnectivityArray_right_copy=0
              Topology_copy=0
              if (molecules_counter.gt.amountofmolecules) then
              else
              !Add 'new first Si atom' to copy matrix/arrays)
              StopCrit3=0
              i=0
              do while (StopCrit3.eq.0)
                  i=i+1
                  if ((ConnectivityMatrix(i,5)+ConnectivityMatrix(i,6)+ConnectivityMatrix(i,7)+ConnectivityMatrix(i,8)).ne.0) then
                      StopCrit3=1
                  endif
              enddo
              do j=1,4*counter_ConnectivityArray,1
                  if (ConnectivityArray_left(j).eq.i) then
                      ConnectivityArray_left_copy(j)=1
                  endif
                  if (ConnectivityArray_right(j).eq.i) then
                      ConnectivityArray_right_copy(j)=1
                  endif
              enddo              
              j=1
              endif
                          else
                              if (counter.ne.1) then  !go back in path array to see if previous CP is still connected to other CPs
                                  path(counter)=0
                                  counter=counter-1
                                  i=path(counter)
                                  j=1
                              else
                                  path(counter)=0 !molecule identified
                                  counter=0
                                  
              if (molecules_counter.eq.(1+(counter_visualization-1)*ceiling(dble(amountofmolecules)/dble(molecules_output)))) then
                                  
                            !      write (5000+molecules_counter,59) Transpose(Configuration)
                                  allocate (Configuration_output(counter_rowConfiguration-1,9))
                                  Configuration_output=0
                                  do i=1,counter_rowConfiguration-1,1
                                      Configuration_output(i,7)=Configuration(i,3)
                                      Configuration_output(i,8)=Configuration(i,4)
                                      Configuration_output(i,9)=Configuration(i,5)
                                  if (i.gt.1) then
                                      include 'SiO.for'
                                      Configuration_output(i,1)=muselected2
                                      Configuration_output(i,2)=1
                                  endif
                                  if (i.gt.2) then
                                  if ((Configuration(Configuration(i,2),1).eq.4).AND.(Configuration(Configuration(i,3),1).eq.3)
     &.AND.(Configuration(Configuration(i,4),1).eq.4)) then
                                      include 'SiOSi.for'
                                      Configuration_output(i,3)=muselected2
                                  else
                                      include 'OSiO.for'
                                      Configuration_output(i,3)=muselected2
                                  endif
                                      Configuration_output(i,4)=1
                                  endif                  
                                  if (i.gt.3) then
                                      include 'SiOSiO.for'
                                      Configuration_output(i,5)=muselected2
                                      Configuration_output(i,6)=1
                                  endif     

                                  enddo
                                  
                                  !Adaption of configuration in case of intramolecular reactions
                                  if (unimolecular_counter.ne.0) then
                                   !   include 'unimolecular.for'
                                  endif
                                  unimolecular_counter=0  !amount of intramolecular paths per molecule
                                  unimolecular_position=0 !position of rows in Unimolecular matrix representing intramolecular paths per molecule
                                  unimolecular_molecule=0 !to make sure we have all intramolecular paths for a certain molecule only once

                                  
                                  write (6000+100*PopulationTimer3+counter_visualization,60) Transpose(Configuration_output)
                                  deallocate (Configuration_output)
                                  allocate (Configuration_type(counter_rowConfiguration-1))
                                  do i=1,counter_rowConfiguration-1,1
                                      if (Configuration(i,1).eq.1) then
                                          Configuration_type(i)='B'
                                      else
                                      if (Configuration(i,1).eq.2) then
                                          Configuration_type(i)='C'
                                      else
                                      if (Configuration(i,1).eq.3) then
                                          Configuration_type(i)='O'
                                      else
                                          Configuration_type(i)='Si'
                                      endif
                                      endif
                                      endif
                                  enddo
                                  write (7000+100*PopulationTimer3+counter_visualization,61) Configuration_type
                                  deallocate (Configuration_type)
                                  
                                  counter_visualization=counter_visualization+1
                                  endif
                                  
                                  molecules_counter=molecules_counter+1
                                  Configuration=0
                                  Configuration(1,1)=4
                                  Configuration(1,2)=1
                                  counter_rowConfiguration=2
              ConnectivityArray_left_copy=0
              ConnectivityArray_right_copy=0
              Topology_copy=0
              if (molecules_counter.gt.amountofmolecules) then
              else
              !Add 'new first Si atom' to copy matrix/arrays)
              StopCrit3=0
              i=0
              do while (StopCrit3.eq.0)
                  i=i+1
                  if ((ConnectivityMatrix(i,5)+ConnectivityMatrix(i,6)+ConnectivityMatrix(i,7)+ConnectivityMatrix(i,8)).ne.0) then
                      StopCrit3=1
                  endif
              enddo
              do j=1,4*counter_ConnectivityArray,1
                  if (ConnectivityArray_left(j).eq.i) then
                      ConnectivityArray_left_copy(j)=1
                  endif
                  if (ConnectivityArray_right(j).eq.i) then
                      ConnectivityArray_right_copy(j)=1
                  endif
              enddo              
              j=1
              endif
                                  
                              endif
                          endif
                      else
                      if (((ConnectivityMatrix(i,1)+ConnectivityMatrix(i,2)+ConnectivityMatrix(i,3)+ConnectivityMatrix(i,4)).eq.0)
     &.AND.((ConnectivityMatrix(i,5)+ConnectivityMatrix(i,6)+ConnectivityMatrix(i,7)+ConnectivityMatrix(i,8)).ne.0)) then
                          do k=1,4,1
                              if (ConnectivityMatrix(i,k+4).ne.0) then
                                  if (counter_rowConfiguration.eq.2) then
                                      Configuration(2,1)=Topology(ConnectivityMatrix(i,k+4),1)
                                      Configuration(2,2)=2
                                      Configuration(2,3)=1
                                      Topology_copy(ConnectivityMatrix(i,k+4),1)=2
                                  else
                                  if (counter_rowConfiguration.eq.3) then
                                      Configuration(3,1)=Topology(ConnectivityMatrix(i,k+4),1)
                                      Configuration(3,2)=3
                                      Configuration(3,3)=1
                                      Configuration(3,4)=2
                                      Topology_copy(ConnectivityMatrix(i,k+4),1)=3
                                  else
                                      Configuration(counter_rowConfiguration,1)=Topology(ConnectivityMatrix(i,k+4),1)
                                      Topology_copy(ConnectivityMatrix(i,k+4),1)=counter_rowConfiguration
                                      Configuration(counter_rowConfiguration,2)=counter_rowConfiguration
                                  Configuration(counter_rowConfiguration,3)=ConnectivityArray_left_copy(ConnectivityMatrix(i,k+4))
                      
                                          if (ConnectivityArray_left_copy(ConnectivityMatrix(i,k+4)).eq.1) then
                                              Configuration(counter_rowConfiguration,4)=2
                                              Configuration(counter_rowConfiguration,5)=3
                                          else
                                          if (ConnectivityArray_left_copy(ConnectivityMatrix(i,k+4)).eq.2) then
                                              Configuration(counter_rowConfiguration,4)=1
                                              Configuration(counter_rowConfiguration,5)=3
                                          else
                                              Configuration(counter_rowConfiguration,4)
     &=Configuration(ConnectivityArray_left_copy(ConnectivityMatrix(i,k+4)),3)
                                              Configuration(counter_rowConfiguration,5)
     &=Configuration(ConnectivityArray_left_copy(ConnectivityMatrix(i,k+4)),4)
                                          endif
                                          endif
                      
                              endif
                              endif

                              counter_rowConfiguration=counter_rowConfiguration+1
                              ConnectivityMatrix(i,k+4)=0
                              endif
                              enddo
                          
                          if (counter.eq.0) then
                          
              if (molecules_counter.eq.(1+(counter_visualization-1)*ceiling(dble(amountofmolecules)/dble(molecules_output)))) then
                                !  write (5000+molecules_counter,59) Transpose(Configuration)
                                  allocate (Configuration_output(counter_rowConfiguration-1,9))
                                  Configuration_output=0
                                  do i=1,counter_rowConfiguration-1,1
                                      Configuration_output(i,7)=Configuration(i,3)
                                      Configuration_output(i,8)=Configuration(i,4)
                                      Configuration_output(i,9)=Configuration(i,5)
                                  if (i.gt.1) then
                                      include 'SiO.for'
                                      Configuration_output(i,1)=muselected2
                                      Configuration_output(i,2)=1
                                  endif
                                  if (i.gt.2) then
                                  if ((Configuration(Configuration(i,2),1).eq.4).AND.(Configuration(Configuration(i,3),1).eq.3)
     &.AND.(Configuration(Configuration(i,4),1).eq.4)) then
                                      include 'SiOSi.for'
                                      Configuration_output(i,3)=muselected2
                                  else
                                      include 'OSiO.for'
                                      Configuration_output(i,3)=muselected2
                                  endif
                                      Configuration_output(i,4)=1
                                  endif                  
                                  if (i.gt.3) then
                                      include 'SiOSiO.for'
                                      Configuration_output(i,5)=muselected2
                                      Configuration_output(i,6)=1
                                  endif     

                                  enddo
                                  
                                  !Adaption of configuration in case of intramolecular reactions
                                  if (unimolecular_counter.ne.0) then
                                    !  include 'unimolecular.for'
                                  endif
                                  unimolecular_counter=0  !amount of intramolecular paths per molecule
                                  unimolecular_position=0 !position of rows in Unimolecular matrix representing intramolecular paths per molecule
                                  unimolecular_molecule=0 !to make sure we have all intramolecular paths for a certain molecule only once

                                  
                                  write (6000+100*PopulationTimer3+counter_visualization,60) Transpose(Configuration_output)
                                  deallocate (Configuration_output)
                                  allocate (Configuration_type(counter_rowConfiguration-1))
                                  do i=1,counter_rowConfiguration-1,1
                                      if (Configuration(i,1).eq.1) then
                                          Configuration_type(i)='B'
                                      else
                                      if (Configuration(i,1).eq.2) then
                                          Configuration_type(i)='C'
                                      else
                                      if (Configuration(i,1).eq.3) then
                                          Configuration_type(i)='O'
                                      else
                                          Configuration_type(i)='Si'
                                      endif
                                      endif
                                      endif
                                  enddo
                                  write (7000+100*PopulationTimer3+counter_visualization,61) Configuration_type
                                  deallocate (Configuration_type)
                                  
                                  counter_visualization=counter_visualization+1
                                  endif
                                  
                                  molecules_counter=molecules_counter+1
                                  Configuration=0
                                  Configuration(1,1)=4
                                  Configuration(1,2)=1
                                  counter_rowConfiguration=2
              ConnectivityArray_left_copy=0
              ConnectivityArray_right_copy=0
              Topology_copy=0
              if (molecules_counter.gt.amountofmolecules) then
              else
              !Add 'new first Si atom' to copy matrix/arrays)
              StopCrit3=0
              i=0
              do while (StopCrit3.eq.0)
                  i=i+1
                  if ((ConnectivityMatrix(i,5)+ConnectivityMatrix(i,6)+ConnectivityMatrix(i,7)+ConnectivityMatrix(i,8)).ne.0) then
                      StopCrit3=1
                  endif
              enddo
              do j=1,4*counter_ConnectivityArray,1
                  if (ConnectivityArray_left(j).eq.i) then
                      ConnectivityArray_left_copy(j)=1
                  endif
                  if (ConnectivityArray_right(j).eq.i) then
                      ConnectivityArray_right_copy(j)=1
                  endif
              enddo              
              j=1
              endif
                          else
                              if (counter.ne.1) then  !go back in path array
                                  path(counter)=0
                                  counter=counter-1
                                  i=path(counter)
                                  j=1
                              else
                                  path(counter)=0
                                  counter=0
                                  
              if (molecules_counter.eq.(1+(counter_visualization-1)*ceiling(dble(amountofmolecules)/dble(molecules_output)))) then
                                  
                              !    write (5000+molecules_counter,59) Transpose(Configuration)
                                  
                                  allocate (Configuration_output(counter_rowConfiguration-1,9))
                                  Configuration_output=0
                                  do i=1,counter_rowConfiguration-1,1
                                      Configuration_output(i,7)=Configuration(i,3)
                                      Configuration_output(i,8)=Configuration(i,4)
                                      Configuration_output(i,9)=Configuration(i,5)
                                  if (i.gt.1) then
                                      include 'SiO.for'
                                      Configuration_output(i,1)=muselected2
                                      Configuration_output(i,2)=1
                                  endif
                                  if (i.gt.2) then
                                  if ((Configuration(Configuration(i,2),1).eq.4).AND.(Configuration(Configuration(i,3),1).eq.3)
     &.AND.(Configuration(Configuration(i,4),1).eq.4)) then
                                      include 'SiOSi.for'
                                      Configuration_output(i,3)=muselected2
                                  else
                                      include 'OSiO.for'
                                      Configuration_output(i,3)=muselected2
                                  endif
                                      Configuration_output(i,4)=1
                                  endif                  
                                  if (i.gt.3) then
                                      include 'SiOSiO.for'
                                      Configuration_output(i,5)=muselected2
                                      Configuration_output(i,6)=1
                                  endif     

                                  enddo
                                  
                                  !Adaption of configuration in case of intramolecular reactions
                                  if (unimolecular_counter.ne.0) then
                                     ! include 'unimolecular.for'
                                  endif
                                  unimolecular_counter=0  !amount of intramolecular paths per molecule
                                  unimolecular_position=0 !position of rows in Unimolecular matrix representing intramolecular paths per molecule
                                  unimolecular_molecule=0 !to make sure we have all intramolecular paths for a certain molecule only once

                                  
                                  write (6000+100*PopulationTimer3+counter_visualization,60) Transpose(Configuration_output)
                                  deallocate (Configuration_output)
                                  allocate (Configuration_type(counter_rowConfiguration-1))
                                  do i=1,counter_rowConfiguration-1,1
                                      if (Configuration(i,1).eq.1) then
                                          Configuration_type(i)='B'
                                      else
                                      if (Configuration(i,1).eq.2) then
                                          Configuration_type(i)='C'
                                      else
                                      if (Configuration(i,1).eq.3) then
                                          Configuration_type(i)='O'
                                      else
                                          Configuration_type(i)='Si'
                                      endif
                                      endif
                                      endif
                                  enddo
                                  write (7000+100*PopulationTimer3+counter_visualization,61) Configuration_type
                                  deallocate (Configuration_type)
                                  
                                  counter_visualization=counter_visualization+1
                                  endif
                                  
                                  molecules_counter=molecules_counter+1
                                  Configuration=0
                                  Configuration(1,1)=4
                                  Configuration(1,2)=1
                                  counter_rowConfiguration=2
              ConnectivityArray_left_copy=0
              ConnectivityArray_right_copy=0
              Topology_copy=0
              if (molecules_counter.gt.amountofmolecules) then
              else
              !Add 'new first Si atom' to copy matrix/arrays)
              StopCrit3=0
              i=0
              do while (StopCrit3.eq.0)
                  i=i+1
                  if ((ConnectivityMatrix(i,5)+ConnectivityMatrix(i,6)+ConnectivityMatrix(i,7)+ConnectivityMatrix(i,8)).ne.0) then
                      StopCrit3=1
                  endif
              enddo
              do j=1,4*counter_ConnectivityArray,1
                  if (ConnectivityArray_left(j).eq.i) then
                      ConnectivityArray_left_copy(j)=1
                  endif
                  if (ConnectivityArray_right(j).eq.i) then
                      ConnectivityArray_right_copy(j)=1
                  endif
              enddo              
              j=1
              endif
                              endif
                          endif  
                      else
                      if (((ConnectivityMatrix(i,1)+ConnectivityMatrix(i,2)+ConnectivityMatrix(i,3)+ConnectivityMatrix(i,4)).ne.0)
     &.AND.((ConnectivityMatrix(i,5)+ConnectivityMatrix(i,6)+ConnectivityMatrix(i,7)+ConnectivityMatrix(i,8)).ne.0)) then
                          do k=1,4,1
                              if ((ConnectivityMatrix(i,k).eq.0).AND.(ConnectivityMatrix(i,k+4).ne.0)) then   !loose segment
                                  
                                  if (counter_rowConfiguration.eq.2) then
                                      Configuration(2,1)=Topology(ConnectivityMatrix(i,k+4),1)
                                      Configuration(2,2)=2
                                      Configuration(2,3)=1
                                      Topology_copy(ConnectivityMatrix(i,k+4),1)=2
                                  else
                                  if (counter_rowConfiguration.eq.3) then
                                      Configuration(3,1)=Topology(ConnectivityMatrix(i,k+4),1)
                                      Configuration(3,2)=3
                                      Configuration(3,3)=1
                                      Configuration(3,4)=2
                                      Topology_copy(ConnectivityMatrix(i,k+4),1)=3
                                  else
                                      Configuration(counter_rowConfiguration,1)=Topology(ConnectivityMatrix(i,k+4),1)
                                      Topology_copy(ConnectivityMatrix(i,k+4),1)=counter_rowConfiguration
                                      Configuration(counter_rowConfiguration,2)=counter_rowConfiguration
                                  Configuration(counter_rowConfiguration,3)=ConnectivityArray_left_copy(ConnectivityMatrix(i,k+4))
                      
                                      if (ConnectivityArray_left_copy(ConnectivityMatrix(i,k+4)).eq.1) then
                                          Configuration(counter_rowConfiguration,4)=2
                                          Configuration(counter_rowConfiguration,5)=3
                                      else
                                      if (ConnectivityArray_left_copy(ConnectivityMatrix(i,k+4)).eq.2) then
                                          Configuration(counter_rowConfiguration,4)=1
                                          Configuration(counter_rowConfiguration,5)=3
                                      else
                                          Configuration(counter_rowConfiguration,4)
     &=Configuration(ConnectivityArray_left_copy(ConnectivityMatrix(i,k+4)),3)
                                  Configuration(counter_rowConfiguration,5)
     &=Configuration(ConnectivityArray_left_copy(ConnectivityMatrix(i,k+4)),4)
                                      endif
                                      endif
                      
                              endif
                              endif

                              counter_rowConfiguration=counter_rowConfiguration+1
                              ConnectivityMatrix(i,k+4)=0                                  
                                  
                              endif
                          enddo
                          
                          l=1
                          do while (ConnectivityMatrix(i,l).eq.0)
                              l=l+1
                          enddo

                          if (counter.eq.0) then
                              path(1)=i
                              path(2)=ConnectivityMatrix(i,l)
                              counter=2
                              
                          else
                              counter=counter+1
                              path(counter)=ConnectivityMatrix(i,l)

                          endif
                          
                          q=1
                          stopcrit3=0
                          do while ((stopcrit3.eq.0).AND.(q.lt.counter))

                              if (path(q).eq.path(counter)) then
                                  stopcrit3=1
                              endif
                                    q=q+1
                          enddo
                          
                          if (stopcrit3.eq.1) then
                          
                              !we have to add O in between
                          
                              Configuration(counter_rowConfiguration,1)=3 !O
                              Configuration(counter_rowConfiguration,2)=counter_rowConfiguration
                              Configuration(counter_rowConfiguration,3)=ConnectivityArray_left_copy(ConnectivityMatrix(i,j+4))
                              Topology_copy(ConnectivityMatrix(i,j+4),1)=counter_rowConfiguration
                              if (Configuration(counter_rowConfiguration,3).eq.1) then
                                  Configuration(counter_rowConfiguration,4)=2
                                  Configuration(counter_rowConfiguration,5)=3
                              else
                              if (Configuration(counter_rowConfiguration,3).eq.2) then
                                  Configuration(counter_rowConfiguration,4)=1
                                  Configuration(counter_rowConfiguration,5)=3
                              else
                                  Configuration(counter_rowConfiguration,4)
     &=Configuration(Configuration(counter_rowConfiguration,3),3)
                                  Configuration(counter_rowConfiguration,5)                                  
     &=Configuration(Configuration(counter_rowConfiguration,3),4)      
                              endif
                              endif
                              counter_rowConfiguration=counter_rowConfiguration+1
                              
                              path(counter)=0
                              counter=counter-1
                          
                          m=1
                          do while (ConnectivityMatrix(ConnectivityMatrix(i,l),m).ne.i) 
                              m=m+1
                          enddo
                          
                          n=i
                          
                          ConnectivityMatrix(ConnectivityMatrix(i,l),m+4)=0
                          ConnectivityMatrix(ConnectivityMatrix(i,l),m)=0

                          i=n
                          j=1
                          
                          ConnectivityMatrix(n,l)=0
                          ConnectivityMatrix(n,l+4)=0
                              
                              
                              
                          else
                              
                          
                  !check if part of loop, go through unimolecular matrix (see distancerule.for) and check whether the new CP is in one of the intramolecular paths
                  x=1
                  do while (Unimolecular(x,1).ne.0)
                  y=1
                      do while (Unimolecular(x,y).ne.0)
                          if ((Unimolecular(x,y).eq.ConnectivityMatrix(i,l)).AND.(Unimolecular_molecule(x).eq.0)) then
                              unimolecular_counter=unimolecular_counter+1
                              unimolecular_position(unimolecular_counter)=x
                              unimolecular_molecule(x)=molecules_counter
                          endif
                          y=y+1
                      enddo
                      x=x+1
                  enddo
                          
                          if (counter_rowConfiguration.eq.2) then
                              Configuration(2,1)=3    !O
                              Configuration(2,2)=2
                              Configuration(2,3)=1
                              Topology_copy(ConnectivityMatrix(i,l+4),1)=2
                              Configuration(3,1)=4    !Si
                              Configuration(3,2)=3
                              Configuration(3,3)=2
                              Configuration(3,4)=1
                              if (ConnectivityArray_left_copy(ConnectivityMatrix(i,l+4)).eq.0) then
                                  ConnectivityArray_left_copy(ConnectivityMatrix(i,l+4))=3
                                  do m=1,4*counter_ConnectivityArray,1
                                      if (ConnectivityArray_left(m).eq.ConnectivityArray_left(ConnectivityMatrix(i,l+4))) then
                                          ConnectivityArray_left_copy(m)=3
                                      else
                                      if (ConnectivityArray_right(m).eq.ConnectivityArray_left(ConnectivityMatrix(i,l+4))) then
                                          ConnectivityArray_right_copy(m)=3
                                      endif
                                      endif
                                  enddo                                  
                             else     
                                  ConnectivityArray_right_copy(ConnectivityMatrix(i,l+4))=3
                                  do m=1,4*counter_ConnectivityArray,1
                                      if (ConnectivityArray_left(m).eq.ConnectivityArray_right(ConnectivityMatrix(i,l+4))) then
                                          ConnectivityArray_left_copy(m)=3
                                      else
                                      if (ConnectivityArray_right(m).eq.ConnectivityArray_right(ConnectivityMatrix(i,l+4))) then
                                          ConnectivityArray_right_copy(m)=3
                                      endif
                                      endif
                                  enddo
                              endif
                              counter_rowConfiguration=4
                          else
                          if (counter_rowConfiguration.eq.3) then
                              Configuration(3,1)=3
                              Configuration(3,2)=3
                              Configuration(3,3)=2
                              Configuration(3,4)=1
                              Topology_copy(ConnectivityMatrix(i,l+4),1)=3
                              Configuration(4,1)=4
                              Configuration(4,2)=4
                              Configuration(4,3)=3
                              Configuration(4,4)=2
                              Configuration(4,5)=1
                              if (ConnectivityArray_left_copy(ConnectivityMatrix(i,l+4)).eq.0) then
                                  ConnectivityArray_left_copy(ConnectivityMatrix(i,l+4))=4
                                  do m=1,4*counter_ConnectivityArray,1
                                      if (ConnectivityArray_left(m).eq.ConnectivityArray_left(ConnectivityMatrix(i,l+4))) then
                                          ConnectivityArray_left_copy(m)=4
                                      else
                                      if (ConnectivityArray_right(m).eq.ConnectivityArray_left(ConnectivityMatrix(i,l+4))) then
                                          ConnectivityArray_right_copy(m)=4
                                      endif
                                      endif
                                  enddo                                  
                             else     
                                  ConnectivityArray_right_copy(ConnectivityMatrix(i,l+4))=4
                                  do m=1,4*counter_ConnectivityArray,1
                                      if (ConnectivityArray_left(m).eq.ConnectivityArray_right(ConnectivityMatrix(i,l+4))) then
                                          ConnectivityArray_left_copy(m)=4
                                      else
                                      if (ConnectivityArray_right(m).eq.ConnectivityArray_right(ConnectivityMatrix(i,l+4))) then
                                          ConnectivityArray_right_copy(m)=4
                                      endif
                                      endif
                                  enddo
                              endif
                              counter_rowConfiguration=5
                          else
                              Configuration(counter_rowConfiguration,1)=3
                              Configuration(counter_rowConfiguration,2)=counter_rowConfiguration
                              if (ConnectivityArray_left_copy(ConnectivityMatrix(i,l+4)).ne.0) then
                                  Configuration(counter_rowConfiguration,3)=ConnectivityArray_left_copy(ConnectivityMatrix(i,l+4))
                              else
                                  Configuration(counter_rowConfiguration,3)=ConnectivityArray_right_copy(ConnectivityMatrix(i,l+4))
                              endif
                              Topology_copy(ConnectivityMatrix(i,l+4),1)=counter_rowConfiguration
                              if (Configuration(counter_rowConfiguration,3).eq.1) then
                                  Configuration(counter_rowConfiguration,4)=2
                                  Configuration(counter_rowConfiguration,5)=3
                              else
                              if (Configuration(counter_rowConfiguration,3).eq.2) then
                                  Configuration(counter_rowConfiguration,4)=1
                                  Configuration(counter_rowConfiguration,5)=3
                              else
                                  Configuration(counter_rowConfiguration,4)
     &=Configuration(Configuration(counter_rowConfiguration,3),3)
                                  Configuration(counter_rowConfiguration,5)                                  
     &=Configuration(Configuration(counter_rowConfiguration,3),4)      
                              endif
                              endif
                              Configuration(counter_rowConfiguration+1,1)=4
                              Configuration(counter_rowConfiguration+1,2)=counter_rowConfiguration+1
                              Configuration(counter_rowConfiguration+1,3)=counter_rowConfiguration
                              Configuration(counter_rowConfiguration+1,4)=Configuration(counter_rowConfiguration,3)
                              Configuration(counter_rowConfiguration+1,5)=Configuration(counter_rowConfiguration,4)
                              
                              if (ConnectivityArray_left_copy(ConnectivityMatrix(i,l+4)).eq.0) then
                                  ConnectivityArray_left_copy(ConnectivityMatrix(i,l+4))=counter_rowConfiguration+1
                                  do m=1,4*counter_ConnectivityArray,1
                                      if (ConnectivityArray_left(m).eq.ConnectivityArray_left(ConnectivityMatrix(i,l+4))) then
                                          ConnectivityArray_left_copy(m)=counter_rowConfiguration+1
                                      else
                                      if (ConnectivityArray_right(m).eq.ConnectivityArray_left(ConnectivityMatrix(i,l+4))) then
                                          ConnectivityArray_right_copy(m)=counter_rowConfiguration+1
                                      endif
                                      endif
                                  enddo                                  
                             else     
                                  ConnectivityArray_right_copy(ConnectivityMatrix(i,l+4))=counter_rowConfiguration+1
                                  do m=1,4*counter_ConnectivityArray,1
                                      if (ConnectivityArray_left(m).eq.ConnectivityArray_right(ConnectivityMatrix(i,l+4))) then
                                          ConnectivityArray_left_copy(m)=counter_rowConfiguration+1
                                      else
                                      if (ConnectivityArray_right(m).eq.ConnectivityArray_right(ConnectivityMatrix(i,l+4))) then
                                          ConnectivityArray_right_copy(m)=counter_rowConfiguration+1
                                      endif
                                      endif
                                  enddo
                              endif
                              counter_rowConfiguration=counter_rowConfiguration+2
                          endif 
                          endif

                          do k=1,4,1
                              if ((ConnectivityMatrix(ConnectivityMatrix(i,l),k).eq.0)
     &.AND.(ConnectivityMatrix(ConnectivityMatrix(i,l),k+4).ne.0)) then   !loose segment
                                  if (counter_rowConfiguration.eq.2) then
                                      Configuration(2,1)=Topology(ConnectivityMatrix(ConnectivityMatrix(i,l),k+4),1)
                                      Configuration(2,2)=2
                                      Configuration(2,3)=1
                                      Topology_copy(ConnectivityMatrix(ConnectivityMatrix(i,l),k+4),1)=2
                                  else
                                  if (counter_rowConfiguration.eq.3) then
                                      Configuration(3,1)=Topology(ConnectivityMatrix(ConnectivityMatrix(i,l),k+4),1)
                                      Configuration(3,2)=3
                                      Configuration(3,3)=1
                                      Configuration(3,4)=2
                                      Topology_copy(ConnectivityMatrix(ConnectivityMatrix(i,l),k+4),1)=3
                                  else
                                      Configuration(counter_rowConfiguration,1)
     &=Topology(ConnectivityMatrix(ConnectivityMatrix(i,l),k+4),1)
                                      Topology_copy(ConnectivityMatrix(ConnectivityMatrix(i,l),k+4),1)=counter_rowConfiguration
                                      Configuration(counter_rowConfiguration,2)=counter_rowConfiguration
                                      Configuration(counter_rowConfiguration,3)
     &=ConnectivityArray_left_copy(ConnectivityMatrix(ConnectivityMatrix(i,l),k+4))
                      
                                      if (ConnectivityArray_left_copy(ConnectivityMatrix(ConnectivityMatrix(i,l),k+4)).eq.1) then
                                          Configuration(counter_rowConfiguration,4)=2
                                          Configuration(counter_rowConfiguration,5)=3
                                      else
                                      if (ConnectivityArray_left_copy(ConnectivityMatrix(ConnectivityMatrix(i,l),k+4)).eq.2) then
                                          Configuration(counter_rowConfiguration,4)=1
                                          Configuration(counter_rowConfiguration,5)=3
                                      else
                                          Configuration(counter_rowConfiguration,4)
     &=Configuration(ConnectivityArray_left_copy(ConnectivityMatrix(ConnectivityMatrix(i,l),k+4)),3)
                                          Configuration(counter_rowConfiguration,5)
     &=Configuration(ConnectivityArray_left_copy(ConnectivityMatrix(ConnectivityMatrix(i,l),k+4)),4)
                                          endif
                                          endif
                                  endif
                                  endif

                                  counter_rowConfiguration=counter_rowConfiguration+1
                                  ConnectivityMatrix(ConnectivityMatrix(i,l),k+4)=0
                              endif
                          enddo
                          
                          m=1
                          do while (ConnectivityMatrix(ConnectivityMatrix(i,l),m).ne.i) 
                              m=m+1
                          enddo
                          
                          n=i
                          
                          ConnectivityMatrix(ConnectivityMatrix(i,l),m+4)=0
                          ConnectivityMatrix(ConnectivityMatrix(i,l),m)=0

                          i=ConnectivityMatrix(i,l)
                          j=1
                          
                          ConnectivityMatrix(n,l)=0
                          ConnectivityMatrix(n,l+4)=0
                  endif        
                  endif
                  endif
                  endif
              endif
      
          else
              if (ConnectivityMatrix(i,j).ne.0) then  !case 3: new CP

                  
                  if (counter.eq.0) then
                      path(1)=i
                      path(2)=ConnectivityMatrix(i,j)
                      counter=2
                  
                  else 
                      counter=counter+1
                      path(counter)=ConnectivityMatrix(i,j)
      
                  endif
                  
                  
                          q=1
                          stopcrit3=0
                          do while ((stopcrit3.eq.0).AND.(q.lt.counter))

                              if (path(q).eq.path(counter)) then
                                  stopcrit3=1
                              endif
                              q=q+1
                          enddo
                          
                          if (stopcrit3.eq.1) then
                          
                              !we have to add O in between
                          
                              Configuration(counter_rowConfiguration,1)=3 !O
                              Configuration(counter_rowConfiguration,2)=counter_rowConfiguration
                              Configuration(counter_rowConfiguration,3)=ConnectivityArray_left_copy(ConnectivityMatrix(i,j+4))
                              Topology_copy(ConnectivityMatrix(i,j+4),1)=counter_rowConfiguration
                              if (Configuration(counter_rowConfiguration,3).eq.1) then
                                  Configuration(counter_rowConfiguration,4)=2
                                  Configuration(counter_rowConfiguration,5)=3
                              else
                              if (Configuration(counter_rowConfiguration,3).eq.2) then
                                  Configuration(counter_rowConfiguration,4)=1
                                  Configuration(counter_rowConfiguration,5)=3
                              else
                                  Configuration(counter_rowConfiguration,4)
     &=Configuration(Configuration(counter_rowConfiguration,3),3)
                                  Configuration(counter_rowConfiguration,5)                                  
     &=Configuration(Configuration(counter_rowConfiguration,3),4)      
                              endif
                              endif
                              counter_rowConfiguration=counter_rowConfiguration+1
                      
                              path(counter)=0
                              counter=counter-1
                              
                              l=1
                              do while (ConnectivityMatrix(ConnectivityMatrix(i,j),l).ne.i) 
                                  l=l+1
                              enddo
                  
                              m=i
                              n=j

                              ConnectivityMatrix(ConnectivityMatrix(i,j),l+4)=0
                              ConnectivityMatrix(ConnectivityMatrix(i,j),l)=0
                  
                              i=m
                              j=1
                  
                              ConnectivityMatrix(m,n)=0
                              ConnectivityMatrix(m,n+4)=0
                                       
                          else
                  
                  !check if part of loop, go through unimolecular matrix (see distancerule.for) and check whether the new CP is in one of the intramolecular paths
                  x=1
                  do while (Unimolecular(x,1).ne.0)
                  y=1
                      do while (Unimolecular(x,y).ne.0)
                          if ((Unimolecular(x,y).eq.ConnectivityMatrix(i,j)).AND.(Unimolecular_molecule(x).eq.0)) then
                              unimolecular_counter=unimolecular_counter+1
                              unimolecular_position(unimolecular_counter)=x
                              unimolecular_molecule(x)=molecules_counter
                          endif
                          y=y+1
                      enddo
                      x=x+1
                  enddo

                  if (counter_rowConfiguration.eq.2) then
                      Configuration(2,1)=3    !O
                      Configuration(2,2)=2
                      Configuration(2,3)=1
                      Topology_copy(ConnectivityMatrix(i,j+4),1)=2
                      Configuration(3,1)=4    !Si
                      Configuration(3,2)=3
                      Configuration(3,3)=2
                      Configuration(3,4)=1
                      if (ConnectivityArray_left_copy(ConnectivityMatrix(i,j+4)).eq.0) then
                          ConnectivityArray_left_copy(ConnectivityMatrix(i,j+4))=3                      
                          do m=1,4*counter_ConnectivityArray,1
                              if (ConnectivityArray_left(m).eq.ConnectivityArray_left(ConnectivityMatrix(i,j+4))) then
                                  ConnectivityArray_left_copy(m)=3
                              else
                              if (ConnectivityArray_right(m).eq.ConnectivityArray_left(ConnectivityMatrix(i,j+4))) then
                                  ConnectivityArray_right_copy(m)=3
                              endif
                              endif
                          enddo
                      else
                          ConnectivityArray_right_copy(ConnectivityMatrix(i,j+4))=3
                          do m=1,4*counter_ConnectivityArray,1
                              if (ConnectivityArray_left(m).eq.ConnectivityArray_right(ConnectivityMatrix(i,j+4))) then
                                  ConnectivityArray_left_copy(m)=3
                              else
                              if (ConnectivityArray_right(m).eq.ConnectivityArray_right(ConnectivityMatrix(i,j+4))) then
                                  ConnectivityArray_right_copy(m)=3
                              endif
                              endif
                          enddo
                      endif                            
                      counter_rowConfiguration=4
                  else
                  if (counter_rowConfiguration.eq.3) then
                      Configuration(3,1)=3
                      Configuration(3,2)=3
                      Configuration(3,3)=2
                      Configuration(3,4)=1
                      Topology_copy(ConnectivityMatrix(i,j+4),1)=3
                      Configuration(4,1)=4
                      Configuration(4,2)=4
                      Configuration(4,3)=3
                      Configuration(4,4)=2
                      Configuration(4,5)=1
                      if (ConnectivityArray_left_copy(ConnectivityMatrix(i,j+4)).eq.0) then
                          ConnectivityArray_left_copy(ConnectivityMatrix(i,j+4))=4
                          do m=1,4*counter_ConnectivityArray,1
                              if (ConnectivityArray_left(m).eq.ConnectivityArray_left(ConnectivityMatrix(i,j+4))) then
                                  ConnectivityArray_left_copy(m)=4
                              else
                              if (ConnectivityArray_right(m).eq.ConnectivityArray_left(ConnectivityMatrix(i,j+4))) then
                                  ConnectivityArray_right_copy(m)=4
                              endif
                              endif
                          enddo                                  
                      else     
                          ConnectivityArray_right_copy(ConnectivityMatrix(i,j+4))=4
                          do m=1,4*counter_ConnectivityArray,1
                              if (ConnectivityArray_left(m).eq.ConnectivityArray_right(ConnectivityMatrix(i,j+4))) then
                                  ConnectivityArray_left_copy(m)=4
                              else
                              if (ConnectivityArray_right(m).eq.ConnectivityArray_right(ConnectivityMatrix(i,j+4))) then
                                  ConnectivityArray_right_copy(m)=4
                              endif
                              endif
                          enddo
                      endif
                      counter_rowConfiguration=5
                  else
                      Configuration(counter_rowConfiguration,1)=3
                      Configuration(counter_rowConfiguration,2)=counter_rowConfiguration
                      if (ConnectivityArray_left_copy(ConnectivityMatrix(i,j+4)).ne.0) then
                          Configuration(counter_rowConfiguration,3)=ConnectivityArray_left_copy(ConnectivityMatrix(i,j+4))
                      else
                          Configuration(counter_rowConfiguration,3)=ConnectivityArray_right_copy(ConnectivityMatrix(i,j+4))
                      endif
                      Topology_copy(ConnectivityMatrix(i,j+4),1)=counter_rowConfiguration
                      if (Configuration(counter_rowConfiguration,3).eq.1) then
                          Configuration(counter_rowConfiguration,4)=2
                          Configuration(counter_rowConfiguration,5)=3
                      else
                      if (Configuration(counter_rowConfiguration,3).eq.2) then
                          Configuration(counter_rowConfiguration,4)=1
                          Configuration(counter_rowConfiguration,5)=3
                      else
                          Configuration(counter_rowConfiguration,4)
     &=Configuration(Configuration(counter_rowConfiguration,3),3)
                          Configuration(counter_rowConfiguration,5)                                  
     &=Configuration(Configuration(counter_rowConfiguration,3),4)      
                      endif
                      endif
                      Configuration(counter_rowConfiguration+1,1)=4
                      Configuration(counter_rowConfiguration+1,2)=counter_rowConfiguration+1
                      Configuration(counter_rowConfiguration+1,3)=counter_rowConfiguration
                      Configuration(counter_rowConfiguration+1,4)=Configuration(counter_rowConfiguration,3)
                      Configuration(counter_rowConfiguration+1,5)=Configuration(counter_rowConfiguration,4)
                      
                      if (ConnectivityArray_left_copy(ConnectivityMatrix(i,j+4)).eq.0) then
                          ConnectivityArray_left_copy(ConnectivityMatrix(i,j+4))=counter_rowConfiguration+1
                          do m=1,4*counter_ConnectivityArray,1
                              if (ConnectivityArray_left(m).eq.ConnectivityArray_left(ConnectivityMatrix(i,j+4))) then
                                  ConnectivityArray_left_copy(m)=counter_rowConfiguration+1
                              else
                              if (ConnectivityArray_right(m).eq.ConnectivityArray_left(ConnectivityMatrix(i,j+4))) then
                                  ConnectivityArray_right_copy(m)=counter_rowConfiguration+1
                              endif
                              endif
                          enddo                                  
                      else     
                          ConnectivityArray_right_copy(ConnectivityMatrix(i,j+4))=counter_rowConfiguration+1
                          do m=1,4*counter_ConnectivityArray,1
                              if (ConnectivityArray_left(m).eq.ConnectivityArray_right(ConnectivityMatrix(i,j+4))) then
                                  ConnectivityArray_left_copy(m)=counter_rowConfiguration+1
                              else
                              if (ConnectivityArray_right(m).eq.ConnectivityArray_right(ConnectivityMatrix(i,j+4))) then
                                  ConnectivityArray_right_copy(m)=counter_rowConfiguration+1
                              endif
                              endif
                          enddo
                      endif
                      counter_rowConfiguration=counter_rowConfiguration+2
                  endif 
                  endif
                  
                  do k=1,4,1
                      if ((ConnectivityMatrix(ConnectivityMatrix(i,j),k).eq.0)
     &.AND.(ConnectivityMatrix(ConnectivityMatrix(i,j),k+4).ne.0)) then   !loose segment
                          if (counter_rowConfiguration.eq.2) then
                              Configuration(2,1)=Topology(ConnectivityMatrix(ConnectivityMatrix(i,j),k+4),1)
                              Configuration(2,2)=2
                              Configuration(2,3)=1
                              Topology_copy(ConnectivityMatrix(ConnectivityMatrix(i,j),k+4),1)=2
                          else
                          if (counter_rowConfiguration.eq.3) then
                              Configuration(3,1)=Topology(ConnectivityMatrix(ConnectivityMatrix(i,j),k+4),1)
                              Configuration(3,2)=3
                              Configuration(3,3)=1
                              Configuration(3,4)=2
                              Topology_copy(ConnectivityMatrix(ConnectivityMatrix(i,j),k+4),1)=3
                          else
                              Configuration(counter_rowConfiguration,1)
     &=Topology(ConnectivityMatrix(ConnectivityMatrix(i,j),k+4),1)
                              Topology_copy(ConnectivityMatrix(ConnectivityMatrix(i,j),k+4),1)=counter_rowConfiguration
                              Configuration(counter_rowConfiguration,2)=counter_rowConfiguration
                              Configuration(counter_rowConfiguration,3)
     &=ConnectivityArray_left_copy(ConnectivityMatrix(ConnectivityMatrix(i,j),k+4))
                      
                              if (ConnectivityArray_left_copy(ConnectivityMatrix(ConnectivityMatrix(i,j),k+4)).eq.1) then
                                  Configuration(counter_rowConfiguration,4)=2
                                  Configuration(counter_rowConfiguration,5)=3
                              else
                              if (ConnectivityArray_left_copy(ConnectivityMatrix(ConnectivityMatrix(i,j),k+4)).eq.2) then
                                  Configuration(counter_rowConfiguration,4)=1
                                  Configuration(counter_rowConfiguration,5)=3
                              else
                                  Configuration(counter_rowConfiguration,4)
     &=Configuration(ConnectivityArray_left_copy(ConnectivityMatrix(ConnectivityMatrix(i,j),k+4)),3)
                                  Configuration(counter_rowConfiguration,5)
     &=Configuration(ConnectivityArray_left_copy(ConnectivityMatrix(ConnectivityMatrix(i,j),k+4)),4)
                              endif
                              endif
                          endif
                          endif

                          counter_rowConfiguration=counter_rowConfiguration+1
                          ConnectivityMatrix(ConnectivityMatrix(i,j),k+4)=0
                      endif
                  enddo
              
                  l=1
                  do while (ConnectivityMatrix(ConnectivityMatrix(i,j),l).ne.i) 
                      l=l+1
                  enddo
                  
                  m=i
                  n=j

                  ConnectivityMatrix(ConnectivityMatrix(i,j),l+4)=0
                  ConnectivityMatrix(ConnectivityMatrix(i,j),l)=0
                  
                  i=ConnectivityMatrix(i,j)
                  j=1
                  
                  ConnectivityMatrix(m,n)=0
                  ConnectivityMatrix(m,n+4)=0
                  
                  endif
     
          else    !if (ConnectivityMatrix(i,j).eq.0).AND.(ConnectivityMatrix(i,j+3).eq.0), nothing really happens (just go one further to the right in connectivitymatrix), except if everything is 0 in connectivitymatrix
              if (j.lt.4) then
                  j=j+1
              else 
                  if (((ConnectivityMatrix(i,1)+ConnectivityMatrix(i,2)+ConnectivityMatrix(i,3)+ConnectivityMatrix(i,4)).eq.0)
     &.AND.((ConnectivityMatrix(i,5)+ConnectivityMatrix(i,6)+ConnectivityMatrix(i,7)+ConnectivityMatrix(i,8)).eq.0)) then
                      if (counter.eq.0) then
                          i=i+1
                          j=1
                      else
                          if (counter.ne.1) then  !go back in path array
                              path(counter)=0
                              counter=counter-1
                              i=path(counter)
                              j=1
                          else
                              path(counter)=0
                              counter=0
                              
              if (molecules_counter.eq.(1+(counter_visualization-1)*ceiling(dble(amountofmolecules)/dble(molecules_output)))) then
                                  
                           !   write (5000+molecules_counter,59) Transpose(Configuration)
                                  allocate (Configuration_output(counter_rowConfiguration-1,9))
                                  Configuration_output=0
                                  do i=1,counter_rowConfiguration-1,1
                                      Configuration_output(i,7)=Configuration(i,3)
                                      Configuration_output(i,8)=Configuration(i,4)
                                      Configuration_output(i,9)=Configuration(i,5)
                                  if (i.gt.1) then
                                      include 'SiO.for'
                                      Configuration_output(i,1)=muselected2
                                      Configuration_output(i,2)=1
                                  endif
                                  if (i.gt.2) then
                                  if ((Configuration(Configuration(i,2),1).eq.4).AND.(Configuration(Configuration(i,3),1).eq.3)
     &.AND.(Configuration(Configuration(i,4),1).eq.4)) then
                                      include 'SiOSi.for'
                                      Configuration_output(i,3)=muselected2
                                  else
                                      include 'OSiO.for'
                                      Configuration_output(i,3)=muselected2
                                  endif
                                      Configuration_output(i,4)=1
                                  endif                  
                                  if (i.gt.3) then
                                      include 'SiOSiO.for'
                                      Configuration_output(i,5)=muselected2
                                      Configuration_output(i,6)=1
                                  endif     
                                  enddo
                                  
                                  !Adaption of configuration in case of intramolecular reactions
                                  if (unimolecular_counter.ne.0) then
                                     ! include 'unimolecular.for'
                                  endif
                                  unimolecular_counter=0  !amount of intramolecular paths per molecule
                                  unimolecular_position=0 !position of rows in Unimolecular matrix representing intramolecular paths per molecule
                                  unimolecular_molecule=0 !to make sure we have all intramolecular paths for a certain molecule only once

                                  
                                  write (6000+100*PopulationTimer3+counter_visualization,60) Transpose(Configuration_output)
                                  deallocate (Configuration_output)
                                  allocate (Configuration_type(counter_rowConfiguration-1))
                                  do i=1,counter_rowConfiguration-1,1
                                      if (Configuration(i,1).eq.1) then
                                          Configuration_type(i)='B'
                                      else
                                      if (Configuration(i,1).eq.2) then
                                          Configuration_type(i)='C'
                                      else
                                      if (Configuration(i,1).eq.3) then
                                          Configuration_type(i)='O'
                                      else
                                          Configuration_type(i)='Si'
                                      endif
                                      endif
                                      endif
                                  enddo
                                  write (7000+100*PopulationTimer3+counter_visualization,61) Configuration_type
                                  deallocate (Configuration_type)  
                                  
                                  counter_visualization=counter_visualization+1
                                  endif
                      
                              molecules_counter=molecules_counter+1
                              Configuration=0
                              Configuration(1,1)=4
                              Configuration(1,2)=1
                              counter_rowConfiguration=2
              ConnectivityArray_left_copy=0
              ConnectivityArray_right_copy=0
              Topology_copy=0
              if (molecules_counter.gt.amountofmolecules) then
              else
              !Add 'new first Si atom' to copy matrix/arrays)
              StopCrit3=0
              i=0
              do while (StopCrit3.eq.0)
                  i=i+1
                  if ((ConnectivityMatrix(i,5)+ConnectivityMatrix(i,6)+ConnectivityMatrix(i,7)+ConnectivityMatrix(i,8)).ne.0) then
                      StopCrit3=1
                  endif
              enddo
              do j=1,4*counter_ConnectivityArray,1
                  if (ConnectivityArray_left(j).eq.i) then
                      ConnectivityArray_left_copy(j)=1
                  endif
                  if (ConnectivityArray_right(j).eq.i) then
                      ConnectivityArray_right_copy(j)=1
                  endif
              enddo              
              j=1
              endif
                                  
                          endif
                      endif
              else
                  if (((ConnectivityMatrix(i,1)+ConnectivityMatrix(i,2)+ConnectivityMatrix(i,3)+ConnectivityMatrix(i,4)).eq.0)
     &.AND.((ConnectivityMatrix(i,5)+ConnectivityMatrix(i,6)+ConnectivityMatrix(i,7)+ConnectivityMatrix(i,8)).ne.0)) then !only loose segments
                      do k=1,4,1
                          if (ConnectivityMatrix(i,k+4).ne.0) then
                              if (counter_rowConfiguration.eq.2) then
                                  Configuration(2,1)=Topology(ConnectivityMatrix(i,k+4),1)
                                  Configuration(2,2)=2
                                  Configuration(2,3)=1
                                  Topology_copy(ConnectivityMatrix(i,k+4),1)=2
                              else
                              if (counter_rowConfiguration.eq.3) then
                                  Configuration(3,1)=Topology(ConnectivityMatrix(i,k+4),1)
                                  Configuration(3,2)=3
                                  Configuration(3,3)=1
                                  Configuration(3,4)=2
                                  Topology_copy(ConnectivityMatrix(i,k+4),1)=3
                              else
                                  Configuration(counter_rowConfiguration,1)=Topology(ConnectivityMatrix(i,k+4),1)
                                  Topology_copy(ConnectivityMatrix(i,k+4),1)=counter_rowConfiguration
                                  Configuration(counter_rowConfiguration,2)=counter_rowConfiguration
                                  Configuration(counter_rowConfiguration,3)=ConnectivityArray_left_copy(ConnectivityMatrix(i,k+4))
                      
                                  if (ConnectivityArray_left_copy(ConnectivityMatrix(i,k+4)).eq.1) then
                                      Configuration(counter_rowConfiguration,4)=2
                                      Configuration(counter_rowConfiguration,5)=3
                                  else
                                  if (ConnectivityArray_left_copy(ConnectivityMatrix(i,k+4)).eq.2) then
                                      Configuration(counter_rowConfiguration,4)=1
                                      Configuration(counter_rowConfiguration,5)=3
                                  else
                                      Configuration(counter_rowConfiguration,4)
     &=Configuration(ConnectivityArray_left_copy(ConnectivityMatrix(i,k+4)),3)
                                      Configuration(counter_rowConfiguration,5)
     &=Configuration(ConnectivityArray_left_copy(ConnectivityMatrix(i,k+4)),4)
                                  endif
                                  endif
                      
                              endif
                              endif

                              counter_rowConfiguration=counter_rowConfiguration+1
                              ConnectivityMatrix(i,k+4)=0
                              endif
                              enddo
                          
                          if (counter.eq.0) then
                              i=1
                              j=1
                          else
                              if (counter.ne.1) then  !go back in path array
                                  path(counter)=0
                                  counter=counter-1
                                  i=path(counter)
                                  j=1
                              else
                                  path(counter)=0
                                  counter=0
                                  
              if (molecules_counter.eq.(1+(counter_visualization-1)*ceiling(dble(amountofmolecules)/dble(molecules_output)))) then
              
                              !    write (5000+molecules_counter,59) Transpose(Configuration)
                                  allocate (Configuration_output(counter_rowConfiguration-1,9))
                                  Configuration_output=0
                                  do i=1,counter_rowConfiguration-1,1
                                      Configuration_output(i,7)=Configuration(i,3)
                                      Configuration_output(i,8)=Configuration(i,4)
                                      Configuration_output(i,9)=Configuration(i,5)
                                  if (i.gt.1) then
                                      include 'SiO.for'
                                      Configuration_output(i,1)=muselected2
                                      Configuration_output(i,2)=1
                                  endif
                                  if (i.gt.2) then
                                      Configuration_output(i,4)=1
                                  endif                  
                                  if (i.gt.3) then
                                      include 'SiOSiO.for'
                                      Configuration_output(i,5)=muselected2
                                      Configuration_output(i,6)=1
                                  endif     
                                  if ((Configuration(Configuration(i,2),1).eq.4).AND.(Configuration(Configuration(i,3),1).eq.3)
     &.AND.(Configuration(Configuration(i,4),1).eq.4)) then
                                      include 'SiOSi.for'
                                      Configuration_output(i,3)=muselected2
                                  else
                                      include 'OSiO.for'
                                      Configuration_output(i,3)=muselected2
                                  endif
                                  enddo
                                  
                                  !Adaption of configuration in case of intramolecular reactions
                                  if (unimolecular_counter.ne.0) then
                                     ! include 'unimolecular.for'
                                  endif
                                  unimolecular_counter=0  !amount of intramolecular paths per molecule
                                  unimolecular_position=0 !position of rows in Unimolecular matrix representing intramolecular paths per molecule
                                  unimolecular_molecule=0 !to make sure we have all intramolecular paths for a certain molecule only once

                                  
                                  write (6000+100*PopulationTimer3+counter_visualization,60) Transpose(Configuration_output)
                                  deallocate (Configuration_output)
                                  allocate (Configuration_type(counter_rowConfiguration-1))
                                  do i=1,counter_rowConfiguration-1,1
                                      if (Configuration(i,1).eq.1) then
                                          Configuration_type(i)='B'
                                      else
                                      if (Configuration(i,1).eq.2) then
                                          Configuration_type(i)='C'
                                      else
                                      if (Configuration(i,1).eq.3) then
                                          Configuration_type(i)='O'
                                      else
                                          Configuration_type(i)='Si'
                                      endif
                                      endif
                                      endif
                                  enddo
                                  write (7000+100*PopulationTimer3+counter_visualization,61) Configuration_type
                                  deallocate (Configuration_type)
                                  
                                  counter_visualization=counter_visualization+1
                                  endif
                                  
                                  molecules_counter=molecules_counter+1
                                  Configuration=0
                                  Configuration(1,1)=4
                                  Configuration(1,2)=1
                                  counter_rowConfiguration=2
              ConnectivityArray_left_copy=0
              ConnectivityArray_right_copy=0
              Topology_copy=0
              if (molecules_counter.gt.amountofmolecules) then
              else
              !Add 'new first Si atom' to copy matrix/arrays)
              StopCrit3=0
              i=0
              do while (StopCrit3.eq.0)
                  i=i+1
                  if ((ConnectivityMatrix(i,5)+ConnectivityMatrix(i,6)+ConnectivityMatrix(i,7)+ConnectivityMatrix(i,8)).ne.0) then
                      StopCrit3=1
                  endif
              enddo
              do j=1,4*counter_ConnectivityArray,1
                  if (ConnectivityArray_left(j).eq.i) then
                      ConnectivityArray_left_copy(j)=1
                  endif
                  if (ConnectivityArray_right(j).eq.i) then
                      ConnectivityArray_right_copy(j)=1
                  endif
              enddo              
              j=1
              endif
                              endif
                          endif
                      else
                      if (((ConnectivityMatrix(i,1)+ConnectivityMatrix(i,2)+ConnectivityMatrix(i,3)+ConnectivityMatrix(i,4)).ne.0)
     &.AND.((ConnectivityMatrix(i,5)+ConnectivityMatrix(i,6)+ConnectivityMatrix(i,7)+ConnectivityMatrix(i,8)).ne.0)) then !at least one new CP
                          do k=1,4,1
                              if ((ConnectivityMatrix(i,k).eq.0).AND.(ConnectivityMatrix(i,k+4).ne.0)) then   !loose segment
                                  
                                  if (counter_rowConfiguration.eq.2) then
                                      Configuration(2,1)=Topology(ConnectivityMatrix(i,k+4),1)
                                      Configuration(2,2)=2
                                      Configuration(2,3)=1
                                      Topology_copy(ConnectivityMatrix(i,k+4),1)=2
                                  else
                                  if (counter_rowConfiguration.eq.3) then
                                      Configuration(3,1)=Topology(ConnectivityMatrix(i,k+4),1)
                                      Configuration(3,2)=3
                                      Configuration(3,3)=1
                                      Configuration(3,4)=2
                                      Topology_copy(ConnectivityMatrix(i,k+4),1)=3
                                  else
                                      Configuration(counter_rowConfiguration,1)=Topology(ConnectivityMatrix(i,k+4),1)
                                      Topology_copy(ConnectivityMatrix(i,k+4),1)=counter_rowConfiguration
                                      Configuration(counter_rowConfiguration,2)=counter_rowConfiguration
                                  Configuration(counter_rowConfiguration,3)=ConnectivityArray_left_copy(ConnectivityMatrix(i,k+4))
                      
                                      if (ConnectivityArray_left_copy(ConnectivityMatrix(i,k+4)).eq.1) then
                                          Configuration(counter_rowConfiguration,4)=2
                                          Configuration(counter_rowConfiguration,5)=3
                                      else
                                      if (ConnectivityArray_left_copy(ConnectivityMatrix(i,k+4)).eq.2) then
                                          Configuration(counter_rowConfiguration,4)=1
                                          Configuration(counter_rowConfiguration,5)=3
                                      else
                                          Configuration(counter_rowConfiguration,4)
     &=Configuration(ConnectivityArray_left_copy(ConnectivityMatrix(i,k+4)),3)
                                  Configuration(counter_rowConfiguration,5)
     &=Configuration(ConnectivityArray_left_copy(ConnectivityMatrix(i,k+4)),4)
                                      endif
                                      endif
                      
                              endif
                              endif

                              counter_rowConfiguration=counter_rowConfiguration+1
                              ConnectivityMatrix(i,k+4)=0                                  
                                  
                              endif
                          enddo
                          
                          !new CP
                          l=1
                          do while (ConnectivityMatrix(i,l).eq.0)
                              l=l+1
                          enddo

                          if (counter.eq.0) then
                              path(1)=i
                              path(2)=ConnectivityMatrix(i,l)
                              counter=2
                              
                          else
                              counter=counter+1
                              path(counter)=ConnectivityMatrix(i,l)

                          endif
                          
              
                          q=1
                          stopcrit3=0
                          do while ((stopcrit3.eq.0).AND.(q.lt.counter))

                              if (path(q).eq.path(counter)) then
                                  stopcrit3=1
                              endif
                              q=q+1
                          enddo
                          
                          if (stopcrit3.eq.1) then
                          
                              !we have to add O in between
                          
                              Configuration(counter_rowConfiguration,1)=3 !O
                              Configuration(counter_rowConfiguration,2)=counter_rowConfiguration
                              Configuration(counter_rowConfiguration,3)=ConnectivityArray_left_copy(ConnectivityMatrix(i,j+4))
                              Topology_copy(ConnectivityMatrix(i,j+4),1)=counter_rowConfiguration
                              if (Configuration(counter_rowConfiguration,3).eq.1) then
                                  Configuration(counter_rowConfiguration,4)=2
                                  Configuration(counter_rowConfiguration,5)=3
                              else
                              if (Configuration(counter_rowConfiguration,3).eq.2) then
                                  Configuration(counter_rowConfiguration,4)=1
                                  Configuration(counter_rowConfiguration,5)=3
                              else
                                  Configuration(counter_rowConfiguration,4)
     &=Configuration(Configuration(counter_rowConfiguration,3),3)
                                  Configuration(counter_rowConfiguration,5)                                  
     &=Configuration(Configuration(counter_rowConfiguration,3),4)      
                              endif
                              endif
                              counter_rowConfiguration=counter_rowConfiguration+1
                          
                              path(counter)=0
                              counter=counter-1
                              
                          m=1
                          do while (ConnectivityMatrix(ConnectivityMatrix(i,l),m).ne.i) 
                              m=m+1
                          enddo
                          
                          n=i
                          
                          ConnectivityMatrix(ConnectivityMatrix(i,l),m+4)=0
                          ConnectivityMatrix(ConnectivityMatrix(i,l),m)=0

                          i=n
                          j=1
                          
                          ConnectivityMatrix(n,l)=0
                          ConnectivityMatrix(n,l+4)=0
                              
                              
                          else
                          
                  !check if part of loop, go through unimolecular matrix (see distancerule.for) and check whether the new CP is in one of the intramolecular paths
                  x=1
                  do while (Unimolecular(x,1).ne.0)
                  y=1
                      do while (Unimolecular(x,y).ne.0)
                          if ((Unimolecular(x,y).eq.ConnectivityMatrix(i,l)).AND.(Unimolecular_molecule(x).eq.0)) then
                              unimolecular_counter=unimolecular_counter+1
                              unimolecular_position(unimolecular_counter)=x
                              unimolecular_molecule(x)=molecules_counter
                          endif
                          y=y+1
                      enddo
                      x=x+1
                  enddo
                          
                          if (counter_rowConfiguration.eq.2) then
                              Configuration(2,1)=3
                              Configuration(2,2)=2
                              Configuration(2,3)=1
                              Topology_copy(ConnectivityMatrix(i,l+4),1)=2
                              Configuration(3,1)=4
                              Configuration(3,2)=3
                              Configuration(3,3)=2
                              Configuration(3,4)=1
                              if (ConnectivityArray_left_copy(ConnectivityMatrix(i,l+4)).eq.0) then
                                  ConnectivityArray_left_copy(ConnectivityMatrix(i,l+4))=3
                                  do m=1,4*counter_ConnectivityArray,1
                                      if (ConnectivityArray_left(m).eq.ConnectivityArray_left(ConnectivityMatrix(i,l+4))) then
                                          ConnectivityArray_left_copy(m)=3
                                      else
                                      if (ConnectivityArray_right(m).eq.ConnectivityArray_left(ConnectivityMatrix(i,l+4))) then
                                          ConnectivityArray_right_copy(m)=3
                                      endif
                                      endif
                                  enddo                                  
                             else     
                                  ConnectivityArray_right_copy(ConnectivityMatrix(i,l+4))=3
                                  do m=1,4*counter_ConnectivityArray,1
                                      if (ConnectivityArray_left(m).eq.ConnectivityArray_right(ConnectivityMatrix(i,l+4))) then
                                          ConnectivityArray_left_copy(m)=3
                                      else
                                      if (ConnectivityArray_right(m).eq.ConnectivityArray_right(ConnectivityMatrix(i,l+4))) then
                                          ConnectivityArray_right_copy(m)=3
                                      endif
                                      endif
                                  enddo
                              endif
                              counter_rowConfiguration=4
                          else
                          if (counter_rowConfiguration.eq.3) then
                              Configuration(3,1)=3
                              Configuration(3,2)=3
                              Configuration(3,3)=2
                              Configuration(3,4)=1
                              Topology_copy(ConnectivityMatrix(i,l+4),1)=3
                              Configuration(4,1)=4
                              Configuration(4,2)=4
                              Configuration(4,3)=3
                              Configuration(4,4)=2
                              Configuration(4,5)=1
                              if (ConnectivityArray_left_copy(ConnectivityMatrix(i,l+4)).eq.0) then
                                  ConnectivityArray_left_copy(ConnectivityMatrix(i,l+4))=4
                                  do m=1,4*counter_ConnectivityArray,1
                                      if (ConnectivityArray_left(m).eq.ConnectivityArray_left(ConnectivityMatrix(i,l+4))) then
                                          ConnectivityArray_left_copy(m)=4
                                      else
                                      if (ConnectivityArray_right(m).eq.ConnectivityArray_left(ConnectivityMatrix(i,l+4))) then
                                          ConnectivityArray_right_copy(m)=4
                                      endif
                                      endif
                                  enddo                                  
                             else     
                                  ConnectivityArray_right_copy(ConnectivityMatrix(i,l+4))=4
                                  do m=1,4*counter_ConnectivityArray,1
                                      if (ConnectivityArray_left(m).eq.ConnectivityArray_right(ConnectivityMatrix(i,l+4))) then
                                          ConnectivityArray_left_copy(m)=4
                                      else
                                      if (ConnectivityArray_right(m).eq.ConnectivityArray_right(ConnectivityMatrix(i,l+4))) then
                                          ConnectivityArray_right_copy(m)=4
                                      endif
                                      endif
                                  enddo
                              endif
                              counter_rowConfiguration=5
                          else
                              Configuration(counter_rowConfiguration,1)=3
                              Configuration(counter_rowConfiguration,2)=counter_rowConfiguration
                              if (ConnectivityArray_left_copy(ConnectivityMatrix(i,l+4)).ne.0) then
                                  Configuration(counter_rowConfiguration,3)=ConnectivityArray_left_copy(ConnectivityMatrix(i,l+4))
                              else
                                  Configuration(counter_rowConfiguration,3)=ConnectivityArray_right_copy(ConnectivityMatrix(i,l+4))
                              endif
                              Topology_copy(ConnectivityMatrix(i,l+4),1)=counter_rowConfiguration
                              if (Configuration(counter_rowConfiguration,3).eq.1) then
                                  Configuration(counter_rowConfiguration,4)=2
                                  Configuration(counter_rowConfiguration,5)=3
                              else
                              if (Configuration(counter_rowConfiguration,3).eq.2) then
                                  Configuration(counter_rowConfiguration,4)=1
                                  Configuration(counter_rowConfiguration,5)=3
                              else
                                  Configuration(counter_rowConfiguration,4)
     &=Configuration(Configuration(counter_rowConfiguration,3),3)
                                  Configuration(counter_rowConfiguration,5)                                  
     &=Configuration(Configuration(counter_rowConfiguration,3),4)      
                              endif
                              endif
                              Configuration(counter_rowConfiguration+1,1)=4
                              Configuration(counter_rowConfiguration+1,2)=counter_rowConfiguration+1
                              Configuration(counter_rowConfiguration+1,3)=counter_rowConfiguration
                              Configuration(counter_rowConfiguration+1,4)=Configuration(counter_rowConfiguration,3)
                              Configuration(counter_rowConfiguration+1,5)=Configuration(counter_rowConfiguration,4)
                              
                              if (ConnectivityArray_left_copy(ConnectivityMatrix(i,l+4)).eq.0) then
                                  ConnectivityArray_left_copy(ConnectivityMatrix(i,l+4))=counter_rowConfiguration+1
                                  do m=1,4*counter_ConnectivityArray,1
                                      if (ConnectivityArray_left(m).eq.ConnectivityArray_left(ConnectivityMatrix(i,l+4))) then
                                          ConnectivityArray_left_copy(m)=counter_rowConfiguration+1
                                      else
                                      if (ConnectivityArray_right(m).eq.ConnectivityArray_left(ConnectivityMatrix(i,l+4))) then
                                          ConnectivityArray_right_copy(m)=counter_rowConfiguration+1
                                      endif
                                      endif
                                  enddo                                  
                             else     
                                  ConnectivityArray_right_copy(ConnectivityMatrix(i,l+4))=counter_rowConfiguration+1
                                  do m=1,4*counter_ConnectivityArray,1
                                      if (ConnectivityArray_left(m).eq.ConnectivityArray_right(ConnectivityMatrix(i,l+4))) then
                                          ConnectivityArray_left_copy(m)=counter_rowConfiguration+1
                                      else
                                      if (ConnectivityArray_right(m).eq.ConnectivityArray_right(ConnectivityMatrix(i,l+4))) then
                                          ConnectivityArray_right_copy(m)=counter_rowConfiguration+1
                                      endif
                                      endif
                                  enddo
                              endif
                              counter_rowConfiguration=counter_rowConfiguration+2
                          endif 
                          endif

                          do k=1,4,1
                              if ((ConnectivityMatrix(ConnectivityMatrix(i,l),k).eq.0)
     &.AND.(ConnectivityMatrix(ConnectivityMatrix(i,l),k+4).ne.0)) then   !loose segments of new CP
                                  if (counter_rowConfiguration.eq.2) then
                                      Configuration(2,1)=Topology(ConnectivityMatrix(ConnectivityMatrix(i,l),k+4),1)
                                      Configuration(2,2)=2
                                      Configuration(2,3)=1
                                      Topology_copy(ConnectivityMatrix(ConnectivityMatrix(i,l),k+4),1)=2
                                  else
                                  if (counter_rowConfiguration.eq.3) then
                                      Configuration(3,1)=Topology(ConnectivityMatrix(ConnectivityMatrix(i,l),k+4),1)
                                      Configuration(3,2)=3
                                      Configuration(3,3)=1
                                      Configuration(3,4)=2
                                      Topology_copy(ConnectivityMatrix(ConnectivityMatrix(i,l),k+4),1)=3
                                  else
                                      Configuration(counter_rowConfiguration,1)
     &=Topology(ConnectivityMatrix(ConnectivityMatrix(i,l),k+4),1)
                                      Topology_copy(ConnectivityMatrix(ConnectivityMatrix(i,l),k+4),1)=counter_rowConfiguration
                                      Configuration(counter_rowConfiguration,2)=counter_rowConfiguration
                                      Configuration(counter_rowConfiguration,3)
     &=ConnectivityArray_left_copy(ConnectivityMatrix(ConnectivityMatrix(i,l),k+4))
                      
                                      if (ConnectivityArray_left_copy(ConnectivityMatrix(ConnectivityMatrix(i,l),k+4)).eq.1) then
                                          Configuration(counter_rowConfiguration,4)=2
                                          Configuration(counter_rowConfiguration,5)=3
                                      else
                                      if (ConnectivityArray_left_copy(ConnectivityMatrix(ConnectivityMatrix(i,l),k+4)).eq.2) then
                                          Configuration(counter_rowConfiguration,4)=1
                                          Configuration(counter_rowConfiguration,5)=3
                                      else
                                          Configuration(counter_rowConfiguration,4)
     &=Configuration(ConnectivityArray_left_copy(ConnectivityMatrix(ConnectivityMatrix(i,l),k+4)),3)
                                          Configuration(counter_rowConfiguration,5)
     &=Configuration(ConnectivityArray_left_copy(ConnectivityMatrix(ConnectivityMatrix(i,l),k+4)),4)
                                          endif
                                          endif
                                  endif
                                  endif

                                  counter_rowConfiguration=counter_rowConfiguration+1
                                  ConnectivityMatrix(ConnectivityMatrix(i,l),k+4)=0
                              endif
                          enddo
                          
                          m=1
                          do while (ConnectivityMatrix(ConnectivityMatrix(i,l),m).ne.i) 
                              m=m+1
                          enddo
                          
                          n=i
                          
                          ConnectivityMatrix(ConnectivityMatrix(i,l),m+4)=0
                          ConnectivityMatrix(ConnectivityMatrix(i,l),m)=0

                          i=ConnectivityMatrix(i,l)
                          j=1
                          
                          ConnectivityMatrix(n,l)=0
                          ConnectivityMatrix(n,l+4)=0
                  endif
                          
                  endif
                  endif
                  endif
              endif
          
      endif
       endif
      endif
           
      enddo
       
      deallocate (ConnectivityMatrix)
      deallocate (path)
      deallocate (Topology_copy)
      deallocate (ConnectivityArray_left_copy)
      deallocate (ConnectivityArray_right_copy)
      deallocate (Unimolecular_molecule)
      deallocate (unimolecular_position)
      deallocate (Configuration)
