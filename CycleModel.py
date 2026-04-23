

"""
=======================================================================================================================
热力循环分析模型由节点、过程、子循环组成
=======================================================================================================================
"""
class CycleModel:
    """
    总循环模型，该循环需要有热源、冷源、中间循环，中间循环需要允许多工质循环的存在，即多个子循环组
    一个子循环组信息需要包括循环组的所有节点，过程以及子循环
    """
    def __init__(self):
        self.


class infoOfNode:
    """
    节点参数，为了方便调用热物性计算模块，该节点参数基本形式等价于MCESP，但存在入流出流等参数，
    """
    def __init__(self, nodeID):
        self.nodeType = None    # 节点等级，0级节点是拓扑节点，计算时不变动
        self.nodeID = nodeID    # 节点编号，便于在其他类中查询编号--int
        self.massFlowIn = 0     # 入流流量--float/kg/s
        self.massFlowOut = 0    # 出流流量--float/kg/s
        self.t = 0              # 温度--float/K
        self.p = 0              # 压力--float/kPa
        self.d = 0              # 密度--float/kg/m^3
        self.h = 0              # 焓--float/kj/kg
        self.s = 0              # 熵--float/kj/kg
        self.cp = 0             # 定压比热容--float/kj/kg/k
        self.ck = 0             # 比热比--float
        self.tc = 0             # 导热系数
        self.vis = 0            # 动力粘度
        self.perOfMedium = None # 工质比例，指定coolprop的单质仅有一个元素的列表--array[]/list[str]
        self.hcantera = None    # 焓值转换备用参数，该参数基于某种工质同样的TP用cantera得到，用于单质的混合--float/kj/kg
        self.RGas = 0           # 气体常数

    def reset(self):
        self.__init__(self.nodeID)

    def setNode(self, propInfo, addInfo=None):
        # 默认输入方法
        for key in propInfo.keys():
            if key in vars(self).keys():
                setattr(self, key, propInfo[key])
        # 附加信息-字典格式
        if addInfo is not None:
            for key in addInfo.keys():
                if key in vars(self).keys():
                    setattr(self, key, addInfo[key])

class infoOfProcess:
    """
    过程参数，该参数包含计算方式（目前只考虑压缩、膨胀和换热）、关联节点、部件效率（叶轮机用等熵效率，换热器用总压恢复系数）、流量
    """
    def __init__(self):
        self.processType = None     # 压缩C、膨胀T、吸热H、放热S
        self.processID = None       # 过程编号

        self.massFlow = 0           # 顺时针过程流量--float/kg/s
        self.efficiency = None      # 等熵效率或总压恢复系数

        self.upstreamNode = None    # 上游节点，记录地址
        self.downstreamNode = None  # 下游节点，记录地址

        self.upstreamNodeID = None    # 上游节点编号
        self.downstreamNodeID = None  # 下游节点编号

        self.cycle = None        # 所在子循环地址
        self.cycleID = None      # 所在子循环编号

class infoOfCycle:
    """
    子循环参数
    """
    def __init__(self):
        self.cycleID = None     # 所在子循环编号
        self.massFlow = None    # 循环流量，用于为各过程流量赋值，可为负数

        self.process = {}       # 过程地址，从左侧顺时针排列，L-U-R-D
        self.processID = {}     # 过程编号

        self.node = {}          # 节点地址，从左下顺时针排列，L-U-R-D
        self.nodeID = {}        # 节点编号

