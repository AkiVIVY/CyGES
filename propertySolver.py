import CoolProp as CP

"+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-"
"物性求解库，仅用coolprop精确求解物质物性"
"+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-"

class propertySolveModel:
    """
    用coolprop精确求解物质物性，导入求解模块中多次调用
    """

    def __init__(self, substanceList):
        """
        初始化时，指定系统计算物质
        :param substanceList--list[str,str, ]
        """
        # -%1%- 创建coolprop物质以及物质类实例列表，后续检测物质类，如果此前调用过，便无需新建，减轻计算负担
        self.substanceList = substanceList  # 记录物质名
        self.substanceClassList = []        # 记录物质接口

        # -%2%- 尝试添加物质，如果不存在，报错
        for substance in self.substanceList:
            try:
                self.substanceClassList.append(CP.AbstractState("HEOS", substance))
            except:
                raise ValueError(
                    '物质名称:' + str(substance) + ';不在Coolprop中，检查propertySolveModel.propertySolve')

    def property(self, inputPar, substance, input1, input2):
        """

        :param input1: 参数1
        :param input2: 参数2
        :param inputPar: 物性求解模式，包括hp, tp, hs, ps--str
        :param substance：工质名--str

        :return: propertyOut：工质物性的字典--dict{str:float}

        coolprop工质名参照
        https://blog.csdn.net/mtc1170824134/article/details/136845572
        http://www.coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids

        常用参数名范围及单位（库的输入输出）：
        字符 单位（CoolProp）  输入输出(Coolprop纯净物)
        T    [K]             IO
        P    [Pa]            IO
        H    [J/kg]          IO
        D    [kg/m^3]        IO
        S    [J/kg/K]        IO
        V    [Pa*s]          O
        C(cp)[J/kg/K]        O
        O(cv)[J/kg/K]        O
        L(tc)[W/m/K]         O
        常用参数对：
        参数对    Coolprop
        HP       gas.update(CoolProp.HmassP_INPUTS, input1*1000, input2*1000)
        TP       gas.update(CoolProp.PT_INPUTS, input2*1000, input1)
        HS       gas.update(CoolProp.HmassSmass_INPUTS, input1*1000, input2*1000)
        PS       gas.update(CoolProp.PSmass_INPUTS, input1*1000, input2*1000)
        目标输出（输入单位一致）：                 输入转换  输出转换
        property['T'] = temperature  # [K]    *1       *1
        property['P'] = pressure  # [kPa]     *1000    /1000
        property['H'] = enthalpy  # [kJ/kg]   *1000    /1000
        property['S'] = entropy  # [kJ/kg/K]  *1000    /1000
        """
        # -%1%- 找到物质索引
        componentNum = self.substanceList.index(substance)

        # -%2%- 根据输入参数，计算物质参数
        if inputPar == 'HP':
            self.substanceClassList[componentNum].update(CP.HmassP_INPUTS, input1 * 1000, input2 * 1000)
        elif inputPar == 'TP':
            self.substanceClassList[componentNum].update(CP.PT_INPUTS, input2 * 1000, input1)
        elif inputPar == 'HS':
            self.substanceClassList[componentNum].update(CP.HmassSmass_INPUTS, input1 * 1000, input2 * 1000)
        elif inputPar == 'PS':
            self.substanceClassList[componentNum].update(CP.PSmass_INPUTS, input1 * 1000, input2 * 1000)
        # -%3%- 创建字典，输出物性
        propertyOut = dict()
        propertyOut['T'] = self.substanceClassList[componentNum].T()
        propertyOut['P'] = self.substanceClassList[componentNum].p() / 1000
        propertyOut['H'] = self.substanceClassList[componentNum].hmass() / 1000
        propertyOut['S'] = self.substanceClassList[componentNum].smass() / 1000
        propertyOut['substance'] = substance
        return propertyOut