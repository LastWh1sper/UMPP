!A unstructured mesh preprocess module
!Mesh reorder
!TODO: mesh division

!Author: BL DONG
!Date: 2020/06/09 

module UnstructuredMeshPreprocess
    
    implicit none
    
    INTEGER :: cellNum,nodeNum
    INTEGER,ALLOCATABLE :: cell(:),cellOld(:)
    INTEGER,ALLOCATABLE :: node(:),nodeOld(:),nodCod(:)
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
INTEGER :: i,j
INTEGER :: nn
INTEGER :: n1,n2,n3
INTEGER :: planPlus
INTEGER :: np1,np2,np3
INTEGER :: oldPlan
OPEN(10, name = "input.mesh")
READ(10,*)nodeNum
ALLOCATE(nodex(nodeNum),nodey(nodeNum),nodCod(nodeNum))
DO i = 1 , nodeNum
    READ(10,*)cache,nodex(i),nodey(i),cache,nodCod(i)
    nodeOld(i) = i
ENDDO

READ(10,*)cellNum
ALLOCATE(cellNode(4,cellNum))

DO i = 1 ,  cellNum
    READ(10,*)cache,(cellNode(j,i),j = 1,4)
    if(cellNode(4,i)==cellNode(3,i))then
        cellNodeNum(i) = 3
    else
        cellNodeNum(i) = 4
    endif
ENDDO
!目前已经读取了一个网格周围的节点,指定第一个网格为第一个平面
!同样的第一个网格周围的节点也在第一个平面上
smallCell = MINVAL(cellOld)
planNum = 1
planNodeNum(1) = 3
planCellNum(1) = 1
cellPlan(smallCell) = 1
do i =1,3

    nn = cellNode(i,planNum)
    planNode(i,1)
    nodeplan(nn) = 1

enddo
!finished 为 cellPlan 中的最小值
finished  = 0

Do while(finished >= 1)

    finished = MINVAL(cellPlan)
    oldPlan = planNum
    planNum = planNum + 1
    !从node开始循环
    Do i = 1, nodeCellNum

        if(nodePlan(i) /= 0)then
            CYCLE
        endif
        n1 = cellNode(1,i)
        n2 = cellNode(2,i)
        n3 = cellNode(3,i)        

        if ((nodePlan(n1) + nodePlan(n2) + nodePlan(n3)) == 0)then
            CYCLE
        endif
        !三角形的某个节点在上一个平面上
        if(n1 == oldPlan .or. n2 == oldPlan .or. n3 == oldPlan)then
            planCellNum(planNum) = planCellNum(planNum) + 1
            !一层包含的位置
            plancell(planCellNum(planNum),planNum) = i
            cellPlan(i) = planNum
            !找到一个网格的三个单元所在的平面
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

close(10)

End subroutine readMeshFile


!readMeshFile子程序计算了每个单元，每个节点所在的平面
!这里对排序后的网格进行重输出

Subroutine outputNewMesh
implicit none 
integer :: i, j 
integer :: nc, nn
integer :: oldn,oldc
integer ::  ic
real :: cache = 0
nc = 0 
nn = 0
!网格重排序
do i = 1, planNum
    do j = 1, planCellNum(i)
        nc = nc + 1
        ic = plancell(j,i)
        cell(nc) = ic
    Enddo
enddo
!节点重排序
do i = 1, planNum

    do j =1, planCellNum(i)

        nn = nn + 1
        in = planNode(j,i)
        node(nn) = in
    enddo
enddo 


OPEN(20,file = "output.mesh")


write(20,*)nodeNum
do i = 1, nodeNum
    oldn = node(i)
    WRITE(20,*)i,nodex(oldn),nodey(oldn),cache,nodCod(oldn)
Enddo
WRITE(20,*)cellNum

do i = 1, cellNum

    oldc = cell(i)
    !这里要输出新的节点
    !200609填坑
    write(20,*)i,(cellNode(j,oldc),j = 1,4)

enddo
End subroutine outputNewMesh


end module UnstructuredMeshPreprocess