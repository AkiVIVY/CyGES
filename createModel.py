import propertySolver


class node:
    """
    任何节点的信息
    """

    def __init__(self, propertyDict):
        self.T = propertyDict['T']
        self.P = propertyDict['P']
        self.H = propertyDict['H']
        self.S = propertyDict['S']
        self.substance = propertyDict['substance']


class workFlow:
    def __init__(self, substanceList, input):
        self.propertySolver = propertySolver.propertySolveModel(substanceList)
        self.substanceList = substanceList

    @staticmethod
    def _normalized_key(value, tolerance):
        """
        按容差把浮点数归一化，用于字典分组键。
        """
        return round(round(value / tolerance) * tolerance, 12)

    @staticmethod
    def _validate_levels(levels, level_name):
        if not isinstance(levels, (list, tuple)) or len(levels) < 2:
            raise ValueError(f'{level_name}必须为至少两个元素的list/tuple')
        for x in levels:
            if not 0 <= x <= 1:
                raise ValueError(f'{level_name}中的分位点必须位于[0,1]区间')

    @staticmethod
    def _build_edges(grouped_nodes, sort_key, edge_type):
        """
        在每组节点内，按指定排序键将相邻节点连接成边。
        """
        edge_list = []
        for node_group in grouped_nodes.values():
            if len(node_group) < 2:
                continue
            sorted_nodes = sorted(node_group, key=sort_key)
            for idx in range(len(sorted_nodes) - 1):
                n1 = sorted_nodes[idx]
                n2 = sorted_nodes[idx + 1]
                edge_list.append({
                    'from_node': n1,
                    'to_node': n2,
                    'type': edge_type,
                    'deltaH': n2.H - n1.H,
                    'deltaT': n2.T - n1.T,
                })
        return edge_list

    def extractSubCycles(self, nodeListInP, tolerance=1e-6):
        """
        从P-S网格节点中提取四边形单元，并打包为字典。
        四边形顶点按顺时针返回：left_bottom, left_top, right_top, right_bottom。
        """
        if tolerance <= 0:
            raise ValueError('tolerance必须大于0')

        p_keys = sorted(nodeListInP.keys())
        node_map = {}
        s_key_set = set()

        # 建立(P,S)->node索引，便于快速查找四个角点
        for p_key, node_group in nodeListInP.items():
            for nodex in node_group:
                s_key = self._normalized_key(nodex.S, tolerance)
                node_map[(p_key, s_key)] = nodex
                s_key_set.add(s_key)

        s_keys = sorted(s_key_set)
        sub_cycle_list = []
        sub_cycle_in_p_band = {}
        sub_cycle_in_s_band = {}
        sub_cycle_id = 0

        for p_idx in range(len(p_keys) - 1):
            p_left = p_keys[p_idx]
            p_right = p_keys[p_idx + 1]
            for s_idx in range(len(s_keys) - 1):
                s_low = s_keys[s_idx]
                s_high = s_keys[s_idx + 1]

                n_lb = node_map.get((p_left, s_low))
                n_lt = node_map.get((p_left, s_high))
                n_rt = node_map.get((p_right, s_high))
                n_rb = node_map.get((p_right, s_low))

                # 仅当四个角点都存在时，认为是一个有效四边形
                if not all([n_lb, n_lt, n_rt, n_rb]):
                    continue

                sub_cycle = {
                    'id': sub_cycle_id,
                    'nodes': {
                        'left_bottom': n_lb,
                        'left_top': n_lt,
                        'right_top': n_rt,
                        'right_bottom': n_rb,
                    },
                    'p_band': (p_left, p_right),
                    's_band': (s_low, s_high),
                }
                sub_cycle_list.append(sub_cycle)

                p_band_key = (p_left, p_right)
                s_band_key = (s_low, s_high)
                sub_cycle_in_p_band.setdefault(p_band_key, []).append(sub_cycle)
                sub_cycle_in_s_band.setdefault(s_band_key, []).append(sub_cycle)
                sub_cycle_id += 1

        return {
            'subCycleList': sub_cycle_list,
            'subCycleInPBand': sub_cycle_in_p_band,
            'subCycleInSBand': sub_cycle_in_s_band,
        }

    "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-"
    "根据温度水平和压力水平生成节点node以及边edge"
    "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-"

    def generateNodesAndEdges(self, TRange, PRange, TLevel, PLevel, substance, tolerance=1e-6):
        """
        根据参数边界、分位点，生成全部节点以及不同压力、温度下的边
        节点按熵、压力分类
        边按压力、熵分类

        :param TRange: 温度上下限-list[min,max]
        :param PRange: 压力上下限-list[min,max]
        :param TLevel: 温度分位点-list[0,x,x,1]
        :param PLevel: 压力分位点-list[0,x,x,1]
        :param substance: 该闭式循环的工质-str
        :param tolerance: 浮点分组/去重容差-float

        :return: dict{
            'nodeList': list[node],
            'nodeListInP': dict[float, list[node]],
            'nodeListInS': dict[float, list[node]],
            'edgeListInP': list[dict],
            'edgeListInS': list[dict]
        }
        """
        # 0.输入参数校验
        if substance not in self.substanceList:
            raise ValueError(f'工质{substance}不在初始化列表中')
        if not isinstance(TRange, (list, tuple)) or len(TRange) != 2:
            raise ValueError('TRange必须为两个元素的list/tuple')
        if not isinstance(PRange, (list, tuple)) or len(PRange) != 2:
            raise ValueError('PRange必须为两个元素的list/tuple')
        Tmin, Tmax = TRange
        Pmin, Pmax = PRange
        if Tmin >= Tmax:
            raise ValueError('TRange最小值必须小于最大值')
        if Pmin >= Pmax:
            raise ValueError('PRange最小值必须小于最大值')
        self._validate_levels(TLevel, 'TLevel')
        self._validate_levels(PLevel, 'PLevel')
        if tolerance <= 0:
            raise ValueError('tolerance必须大于0')

        # 0.参数准备
        # 参数范围和总差值
        dP = Pmax - Pmin
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
                except Exception as exc:
                    raise ValueError(
                        f'generateNodesAndEdges时，TP参数对超边界: T={Tx}, P={Px}, substance={substance}'
                    ) from exc
        # 2.生成二级节点
        # 记录最低温度和最高温度的熵列表，每跳出一个范围，减少一个二级节点
        SLowList = [self.propertySolver.property('TP', substance, Tmin, Px)['S'] for Px in Plist]
        SHighList = [self.propertySolver.property('TP', substance, Tmax, Px)['S'] for Px in Plist]
        nodeListB = []

        # 对于所有节点，在两个熵范围间的添加所有不相等的压力，在高熵或低熵范围的，削减对应的压力水平
        for nodex in nodeListA:
            for idx, Px in enumerate(Plist):
                # 在不超温的界限内找对应的等熵二级节点
                if not (SLowList[idx] - tolerance <= nodex.S <= SHighList[idx] + tolerance):
                    continue
                if abs(nodex.P - Px) <= tolerance:
                    continue
                try:
                    nodeListB.append(node(self.propertySolver.property('PS', substance, Px, nodex.S)))
                except Exception as exc:
                    raise ValueError(
                        f'generateNodesAndEdges时，PS参数对超边界: S={nodex.S}, P={Px}, substance={substance}'
                    ) from exc

        # 3.将节点分别按压力、熵分离排列，然后排序方式为从小到大，
        nodeListRaw = nodeListA + nodeListB

        # 3.1 基于P/S容差去重
        nodeList = []
        seen = set()
        for nodex in nodeListRaw:
            p_key = self._normalized_key(nodex.P, tolerance)
            s_key = self._normalized_key(nodex.S, tolerance)
            key = (p_key, s_key)
            if key in seen:
                continue
            seen.add(key)
            nodeList.append(nodex)

        # 3.2 按P和S分组
        nodeListInP = {}
        nodeListInS = {}
        for nodex in nodeList:
            p_key = self._normalized_key(nodex.P, tolerance)
            s_key = self._normalized_key(nodex.S, tolerance)
            if p_key in nodeListInP:
                nodeListInP[p_key].append(nodex)
            else:
                nodeListInP[p_key] = [nodex]
            if s_key in nodeListInS:
                nodeListInS[s_key].append(nodex)
            else:
                nodeListInS[s_key] = [nodex]

        # 4.在每个压力和每个熵水平下，将所有相邻两点组成一组边
        edgeListInP = self._build_edges(nodeListInP, sort_key=lambda n: n.S, edge_type='isobaric')
        edgeListInS = self._build_edges(nodeListInS, sort_key=lambda n: n.P, edge_type='isentropic')

        subCycleDict = self.extractSubCycles(nodeListInP, tolerance=tolerance)

        return {
            'nodeList': nodeList,
            'nodeListInP': nodeListInP,
            'nodeListInS': nodeListInS,
            'edgeListInP': edgeListInP,
            'edgeListInS': edgeListInS,
            'subCycleDict': subCycleDict,
        }