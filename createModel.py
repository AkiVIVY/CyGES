import propertySolver

class node:
    """
    任何节点的信息
    """
    def __init__(self, propertyDict):
        self.T = propertyDict['T']
        self.P = propertyDict['P']
        self.H = propertyDict['H']
        self.S = propertyDict['PS']
        self.substance = propertyDict['substance']

class workFlow:
    def __init__(self, substanceList, input):
        self.propertySolver = propertySolver.propertySolveModel(substanceList)


    "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-"
    "根据温度水平和压力水平，生成节点node，边edge"
    "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-"

    def generateNodesAndEdges(self, TRange, PRange, TLevel, PLevel, substance):
        """
        根据参数边界、分位点，生成全部节点以及不同压力、温度下的边
        节点按熵、压力分类
        边按压力、熵分类

        :param TRange: 温度上下限-list[min,max]
        :param PRange: 压力上下限-list[min,max]
        :param TLevel: 温度分位点-list[0,x,x,1]
        :param PLevel: 压力分位点-list[0,x,x,1]
        :param substance: 该闭式循环的工质-str
        """
        # 0.参数准备
        # 参数范围和总差值
        Pmin = PRange[0]
        Pmax = PRange[1]
        dP = Pmax - Pmin
        Tmin = TRange[0]
        Tmax = TRange[1]
        dT = Tmax - Tmin
        # 压力及温度水平列表
        Plist = []
        Tlist = []
        # 根据百分位确定压力及温度水平
        for Px in PLevel:
            Plist.append(Pmin + dP * Px)
        for Tx in TLevel:
            Tlist.append(Tmin + dT * Tx)
        # 1.生成一级节点
        nodeListA = []     # 所有由TP生成的一级节点
        for Px in Plist:
            for Tx in Tlist:
                try:
                    nodeListA.append(node(self.propertySolver.property('TP', substance, Tx, Px)))
                except:
                    raise ValueError('generateNodesAndEdges时，TP参数对超边界')
        # 2.生成二级节点
        # 记录最低温度和最高温度的熵列表，每跳出一个范围，减少一个二级节点
        SLowList = [self.propertySolver.property('TP', substance, Tmin, Px)['S'] for Px in Plist]
        SHighList = [self.propertySolver.property('TP', substance, Tmax, Px)['S'] for Px in Plist]
        nodeListB = []
        # 对于所有节点，在两个熵范围间的添加所有不相等的压力，在高熵或低熵范围的，削减对应的压力水平
        for nodex in nodeListA:
            SLowNum = 0
            SHighNum = 0
            for Sx in SLowList:
                # 仅当节点低于该低温熵，跳过并计数
                if nodex.S < Sx:
                    SLowNum += 1
            for Sx in SHighList:
                # 如果节点熵高于该高温熵，跳过并计数
                if nodex.S > Sx:
                    SHighNum += 1
            # 在不超温的界限内找对应的等熵二级节点
            for Px in Plist[SLowNum: len(Plist) - SHighNum - 1]:
                if nodex.P == Px:
                    pass
                else:
                    nodeListB.append(node(self.propertySolver.property('PS', substance, nodex.S, Px)))
        # 3.将节点分别按压力、熵分离，排序方式为从小到大，
        nodeList = nodeListA + nodeListB
        nodeListInP = {}
        nodeListInS = {}
        for nodex in nodeList:
            if nodex.P in nodeListInP:
                nodeListInP[nodex.P].append(nodex)
            else:
                nodeListInP[nodex.P] = [nodex]
            if nodex.S in nodeListInS:
                nodeListInP[nodex.P].append(nodex)


            nodeListInS[nodex.S].append(nodex)