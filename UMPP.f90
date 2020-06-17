!A unstructured mesh preprocess module
!Mesh reorder
!TODO: mesh division

!Author: BL DONG
!Date: 2020/06/09 

module UnstructuredMeshPreprocess
    
    implicit none
    
    INTEGER :: cellNum,nodeNum
    !cellSta  = 1 分配成功
    INTEGER,ALLOCATABLE :: cell(:),cellOld(:),cellNew(:),cellSta(:)
    INTEGER,ALLOCATABLE :: node(:),nodeOld(:),nodCod(:),nodeNew(:)
    REAL,ALLOCATABLE :: nodex(:),nodey(:)
    INTEGER :: smallCell
    !一个网格周围的网格
    INTEGER,ALLOCATABLE :: nearCell(:,:),nearCellNUm(:)
    !一个网格的节点
    INTEGER,ALLOCATABLE :: cellNode(:,:),cellNodeNum(:)
    !一个节点周围的网格
    INTEGER,ALLOCATABLE :: nodeCell(:,:),nodeCellNum(:)
    !每个cell所在超平面, 每个超平面包含的个数
    INTEGER,ALLOCATABLE :: cellPlan(:),planCell(:,:),planCellNum(:)
    INTEGER,ALLOCATABLE :: nodePlan(:),planNode(:,:),planNodeNum(:)
    !平面数
    INTEGER :: planNum 
    ! 是否分配完成
    INTEGER :: finished

contains

!读取mesh文件
subroutine readMeshFile
implicit none 
REAL :: cache
integer :: cache1
INTEGER :: i,j
INTEGER :: nn
INTEGER :: n1,n2,n3
INTEGER :: planPlus
INTEGER :: np1,np2,np3
INTEGER :: oldPlan 
INTEGER :: countCells = 0
OPEN(10, name = "input.mesh")
READ(10,*)nodeNum
write(*,*)"The mesh file has ", nodeNum, " nodes"
!声明plannode 与 plannode num 最多允许有10000 个层 每个层最多有10000个变量
ALLOCATE(nodex(nodeNum),nodey(nodeNum),nodCod(nodeNum),node(nodeNum),nodeOld(nodeNum),nodeNew(nodeNum),nodePlan(nodeNum),planNode(10000,10000),planNodeNum(10000))
DO i = 1 , nodeNum
    READ(10,*)nodeold(i),nodex(i),nodey(i),cache1,nodCod(i)
ENDDO

Write(*,*)" read nodes"
READ(10,*)cellNum

ALLOCATE(cellNode(4,cellNum),cellNodenum(cellNum),cellOld(cellNum),cell(cellNum),cellNew(cellNum),cellPlan(cellNum),planCell(10000,10000),planCellNum(10000),cellSta(cellNum))
cellSta = 0
nodeplan = 0
cellPlan = 0 
DO i = 1 ,  cellNum
    READ(10,*)cellold(i),(cellNode(j,i),j = 1,4)
    if(cellNode(4,i)==cellNode(3,i))then
        cellNodeNum(i) = 3
    else
        cellNodeNum(i) = 4
    endif
ENDDO

write(*,*)" read mesh" 
!目前已经读取了一个网格周围的节点,指定第一个网格为第一个平面
!同样的第一个网格周围的节点也在第一个平面上
smallCell = MINVAL(cellOld)
planNum = 1
plancellnum = 0
plannodenum = 0
planNodeNum(1) = 3
planCellNum(1) = 1
countCells = 1
cellPlan(smallCell) = 1
cellSta(smallCell)  = 1
plancell(1,1) = smallCell
do i =1,3

    nn = cellNode(i,smallCell)
    planNode(i,1)  = nn
    nodeplan(nn) = 1

enddo
!finished 为 cellPlan 中的最小值
finished  = 0

Do while(finished == 0)

    finished = MINVAL(cellPlan)
    oldPlan = planNum
    planNum = planNum + 1
    
    write(*,*)"plan",oldplan
    write(*,*)"plancellnum",planCellNUm(oldplan)
    !从node开始循环
    Do i = 1, cellNum
        
        
        if(cellPlan(i) /= 0)then
            CYCLE
        endif
        n1 = cellNode(1,i)
        n2 = cellNode(2,i)
        n3 = cellNode(3,i)        

        if ((nodePlan(n1) + nodePlan(n2) + nodePlan(n3)) == 0)then
            CYCLE
        endif
        !三角形的某个节点在上一个平面上
        if(nodeplan(n1) == oldPlan .or. nodeplan(n2) == oldPlan .or. nodeplan(n3) == oldPlan)then
            cellsta(i) = 1
            planCellNum(planNum) = planCellNum(planNum) + 1
            countCells = countCells + 1
            !一层包含的位置
            plancell(planCellNum(planNum),planNum) = i
            cellPlan(i) = planNum
            write(*,*) i,plannum
            !找到一个网格的三个节点所在的平面
            if(nodePlan(n1)==0)then
                planNodeNum(planNum) = planNodeNum(planNum) + 1
                planNode(planNodeNum(planNum),planNum) = n1
                nodePlan(n1) = planNum
            endif

            if(nodePlan(n2)==0)then
                planNodeNum(planNum) = planNodeNum(planNum) + 1
                planNode(planNodeNum(planNum),planNum) = n2
                nodePlan(n2) = planNum
            endif

            if(nodePlan(n3)==0)then
                planNodeNum(planNum) = planNodeNum(planNum) + 1
                planNode(planNodeNum(planNum),planNum) = n3
                nodePlan(n3) = planNum
            endif 
        endif
    Enddo

    finished = minval(cellPlan)

Enddo


Do i  = 1, cellNum
    
    if(cellSta(i)==0)then
        
        n1 = cellNode(1,i)
        n2 = cellNode(2,i)
        n3 = cellNode(3,i)       
        write(*,*)"cell & plan",i
        write(*,*)"nodesplan",nodeplan(n1),nodePlan(n2),nodePlan(n3)
        write(*,*)"cell",cellPlan(i)
    
    endif
    
Enddo


if(countCells /= cellNum)then
    write(*,*)"c0untcells",countcells,cellnum
    write(*,*)"error in finding cells"
endif

close(10)
End subroutine readMeshFile



!readMeshFile子程序计算了每个单元，每个节点所在的平面
!这里对排序后的网格进行重输出

Subroutine outputNewMesh
implicit none 
integer :: i, j 
integer :: nc, nn
integer :: oldn,oldc
integer ::  ic,in
integer :: nnd
integer :: oldnod(4)
integer :: newnod(4)
real :: cache = 0
nc = 0 
nn = 0
!网格重排序
do i = 1, planNum
    do j = 1, planCellNum(i)
        nc = nc + 1
        ic = plancell(j,i)
        cell(nc) = ic
        cellnew(ic) = nc
    Enddo
enddo


!节点重排序
do i = 1, planNum

    do j =1, planNodeNum(i)

        nn = nn + 1
        in = planNode(j,i)

        node(nn) = in
        nodeNew(in) = nn
    enddo
enddo 


OPEN(20,file = "output.mesh")


write(20,*)nodeNum

do i = 1, nodeNum
    oldn = node(i)
    write(*,*)"oldn",oldn
    WRITE(20,*)i,nodex(oldn),nodey(oldn),cache,nodCod(oldn)
Enddo

WRITE(20,*)cellNum

do i = 1, cellNum

    oldc = cell(i)
    
    oldnod(:) = cellNode(:,oldc)
    oldnod(4) = oldnod(3)
    write(*,*)"outputcell",i,oldc,oldnod(:)
    do j =1, 4

        !重排序后的三角形网格对应新的节点
        newnod(j) = nodenew(oldnod(j))

    enddo

    !这里要输出新的节点
    !200609填坑
    write(20,*)i,newnod(:)

enddo
End subroutine outputNewMesh


end module UnstructuredMeshPreprocess