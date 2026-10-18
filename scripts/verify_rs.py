# -*- coding: utf-8 -*-
"""
验证RS构型学习工具的科学正确性
选取315化学(农)真题风格的三道题：基础/中档/难题
"""
import math

# ============ CIP优先级数据（与HTML中CIP_CHILDREN一致）============
ATOMIC_NUM = {'H':1,'D':2,'C':6,'N':7,'O':8,'F':9,'P':15,'S':16,'Cl':17,'Br':35,'I':53}

CIP_CHILDREN = {
    'H':[], 'D':[], 'F':[], 'Cl':[], 'Br':[], 'I':[],
    'CH3':['H','H','H'],
    'C2H5':['CH3','H','H'],
    'nPr':['C2H5','H','H'],
    'iPr':['CH3','CH3','H'],
    'tBu':['CH3','CH3','CH3'],
    'iBu':['iPr','H','H'],
    'sBu':['CH3','C2H5','H'],
    'Vinyl':['=CH2','H'],
    'Allyl':['Vinyl','H','H'],
    'Ethynyl':['≡CH'],
    'Ph':['PhC','PhC','H'],
    'PhC':['PhC','PhCH','H'],
    'PhCH':['PhC','H'],
    '=CH2':['H','H'],
    '≡CH':['H'],
    'Bn':['Ph','H','H'],
    'OH':['H'],
    'OCH3':['CH3'],
    'OC2H5':['C2H5'],
    'OPh':['Ph'],
    'OAc':['COCH3'],
    'CHO':['O','O','H'],
    'COCH3':['O','O','CH3'],
    'COOH':['O','O','OH'],
    'COOCH3':['O','O','OCH3'],
    'CONH2':['O','O','NH2'],
    'CH2OH':['OH','H','H'],
    'CN':['N','N','N'],
    'NH2':['H','H'],
    'NHCH3':['CH3','H'],
    'NMe2':['CH3','CH3'],
    'NO2':['O','O','O'],
    'NHAc':['COCH3','H'],
    'SH':['H'],
    'SCH3':['CH3'],
    'SO3H':['O','O','OH'],
}

GROUP_ATOM = {
    'H':'H','D':'H','F':'F','Cl':'Cl','Br':'Br','I':'I',
    'CH3':'C','C2H5':'C','nPr':'C','iPr':'C','tBu':'C','iBu':'C','sBu':'C',
    'Vinyl':'C','Allyl':'C','Ethynyl':'C','Ph':'C','Bn':'C',
    'OH':'O','OCH3':'O','OC2H5':'O','OPh':'O','OAc':'O',
    'CHO':'C','COCH3':'C','COOH':'C','COOCH3':'C','CONH2':'C','CH2OH':'C','CN':'C',
    'NH2':'N','NHCH3':'N','NMe2':'N','NO2':'N','NHAc':'N',
    'SH':'S','SCH3':'S','SO3H':'S',
    '=CH2':'C','≡CH':'C','PhC':'C','PhCH':'C',
}

cip_cache = {}
def build_cip_layers(key, depth=0, visited=None):
    if depth > 8: return []
    cache_key = key + '_' + str(depth)
    if cache_key in cip_cache: return cip_cache[cache_key]
    if visited is None: visited = set()
    if key in visited: return []
    visited.add(key)
    children = CIP_CHILDREN.get(key, [])
    if not children:
        cip_cache[cache_key] = []
        return []
    layer = sorted([ATOMIC_NUM.get(GROUP_ATOM.get(c, c), 6) for c in children], reverse=True)
    next_layers = []
    for c in children:
        sub = build_cip_layers(c, depth+1, set(visited))
        for i, s in enumerate(sub):
            if i >= len(next_layers): next_layers.append([])
            next_layers[i].extend(s)
    result = [layer] + [sorted(l, reverse=True) for l in next_layers]
    cip_cache[cache_key] = result
    return result

def compare_cip(a, b):
    an = ATOMIC_NUM.get(GROUP_ATOM.get(a, a), 0)
    bn = ATOMIC_NUM.get(GROUP_ATOM.get(b, b), 0)
    if an != bn: return an - bn
    la, lb = build_cip_layers(a), build_cip_layers(b)
    for d in range(max(len(la), len(lb))):
        ra = la[d] if d < len(la) else []
        rb = lb[d] if d < len(lb) else []
        for i in range(max(len(ra), len(rb))):
            va = ra[i] if i < len(ra) else 0
            vb = rb[i] if i < len(rb) else 0
            if va != vb: return va - vb
    return 0

# ============ R/S计算（与HTML中calcRS一致）============
def normalize(v):
    l = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    return [v[0]/l, v[1]/l, v[2]/l] if l > 0 else v

def rotate_point(v, axis, angle):
    x, y, z = v
    ux, uy, uz = axis
    cos_a = math.cos(angle)
    sin_a = math.sin(angle)
    rx = x*cos_a + (uy*z - uz*y)*sin_a + ux*(ux*x+uy*y+uz*z)*(1-cos_a)
    ry = y*cos_a + (uz*x - ux*z)*sin_a + uy*(ux*x+uy*y+uz*z)*(1-cos_a)
    rz = z*cos_a + (ux*y - uy*x)*sin_a + uz*(ux*x+uy*y+uz*z)*(1-cos_a)
    return [rx, ry, rz]

def get_tetra_coords(groups, wedge):
    # 与HTML中一致的真正四面体坐标
    base = [[0,1.25,-0.45],[1.15,0,0.35],[0,-1.25,-0.45],[-1.15,0,0.35]]
    coords = []
    for i in range(4):
        x, y, z = base[i]
        if wedge[i] == 'up': z += 0.85
        elif wedge[i] == 'down': z -= 0.85
        coords.append([x, y, z])
    return coords

def calc_rs(groups, wedge):
    # CIP排序
    idx = [0,1,2,3]
    idx.sort(key=lambda i: -compare_cip(groups[i], groups[0]) if False else 0)
    idx.sort(key=lambda i: compare_cip(groups[i], groups[0]))
    # 正确排序：从高到低
    idx = sorted(range(4), key=lambda i: compare_cip(groups[i], groups[0]), reverse=True)
    # 重新正确排序
    pairs = [(compare_cip(groups[i], groups[0]), i) for i in range(4)]
    # 用两两比较排序
    idx = sorted(range(4), key=lambda i: sum(1 for j in range(4) if compare_cip(groups[i], groups[j]) > 0), reverse=True)
    
    coords = get_tetra_coords(groups, wedge)
    ordered = [coords[i] for i in idx]
    v4 = ordered[3]
    n = normalize(v4)
    target = [0, 0, -1]
    axis = [n[1]*target[2]-n[2]*target[1], n[2]*target[0]-n[0]*target[2], n[0]*target[1]-n[1]*target[0]]
    al = math.sqrt(axis[0]**2 + axis[1]**2 + axis[2]**2)
    angle = math.acos(max(-1, min(1, n[0]*target[0]+n[1]*target[1]+n[2]*target[2])))
    if al < 1e-6:
        r1, r2, r3 = ordered[0], ordered[1], ordered[2]
    else:
        axis = [axis[0]/al, axis[1]/al, axis[2]/al]
        r1 = rotate_point(ordered[0], axis, angle)
        r2 = rotate_point(ordered[1], axis, angle)
        r3 = rotate_point(ordered[2], axis, angle)
    cross2d = (r2[0]-r1[0])*(r3[1]-r1[1]) - (r2[1]-r1[1])*(r3[0]-r1[0])
    return 'S' if cross2d > 0 else 'R', idx

# ============ 验证三道题 ============
print("=" * 60)
print("315化学(农) RS构型真题验证")
print("=" * 60)

# ---- 基础题：2013年第47题 (2R)-2-氯丁醛 ----
print("\n【基础题】2013年第47题：(2R)-2-氯丁醛")
print("结构：OHC-CH(Cl)-CH2-CH3")
print("手性碳C2连接：CHO, Cl, C2H5, H")
groups1 = ['CHO', 'Cl', 'C2H5', 'H']
wedge1 = ['plain', 'up', 'plain', 'down']  # CHO上(纸面), Cl右(实楔前), C2H5下(纸面), H左(虚楔后)
rs1, order1 = calc_rs(groups1, wedge1)
cip_order1 = [groups1[i] for i in order1]
print(f"  CIP优先级：{cip_order1[0]} > {cip_order1[1]} > {cip_order1[2]} > {cip_order1[3]}")
print(f"  工具计算结果：{rs1} 构型")
print(f"  真题答案：R 构型")
print(f"  验证：{'✓ 正确' if rs1 == 'R' else '✗ 错误'}")

# ---- 中档题：Fischer投影式判断（2013年第21题类型）----
print("\n【中档题】Fischer投影式判断（2013年第21题类型）")
print("Fischer投影式：COOH在上, OH在右, CH3在下, H在左")
print("（横前竖后，H在横键朝前）")
# Fischer转楔形键：横键朝前(up)，竖键朝后(down)
groups2 = ['COOH', 'OH', 'CH3', 'H']
wedge2 = ['down', 'up', 'down', 'up']  # 上COOH朝后, 右OH朝前, 下CH3朝后, 左H朝前
rs2, order2 = calc_rs(groups2, wedge2)
cip_order2 = [groups2[i] for i in order2]
print(f"  CIP优先级：{cip_order2[0]} > {cip_order2[1]} > {cip_order2[2]} > {cip_order2[3]}")
print(f"  工具计算结果：{rs2} 构型")
print(f"  手动验证：OH>COOH>CH3>H, H在横键(前), OH→COOH→CH3为逆时针, H在前故逆R → R")
print(f"  正确答案：R 构型")
print(f"  验证：{'✓ 正确' if rs2 == 'R' else '✗ 错误'}")

# ---- 难题：双手性中心 (2R,3R)-酒石酸（2014年第22题类型）----
print("\n【难题】双手性中心判断（2014年第22题类型）")
print("(2R,3R)-酒石酸：HOOC-CH(OH)-CH(OH)-COOH")
print("C2连接：COOH, OH, C3片段, H")
print("C3连接：C2片段, OH, COOH, H")

# C2: top=COOH, right=OH(up), bottom=C3片段, left=H(down)
groups_c2 = ['COOH', 'OH', '_frag_c3', 'H']
wedge_c2 = ['plain', 'up', 'plain', 'down']
# C3片段：连OH, COOH, H → 锚点C连O,C,H
CIP_CHILDREN['_frag_c3'] = ['OH', 'COOH', 'H']
GROUP_ATOM['_frag_c3'] = 'C'
rs_c2, order_c2 = calc_rs(groups_c2, wedge_c2)

# C3: top=C2片段, right=H(down), bottom=COOH, left=OH(up)
# (2R,3R)-酒石酸中C2和C3的OH在相反侧
groups_c3 = ['_frag_c2', 'H', 'COOH', 'OH']
wedge_c3 = ['plain', 'down', 'plain', 'up']
CIP_CHILDREN['_frag_c2'] = ['OH', 'COOH', 'H']
GROUP_ATOM['_frag_c2'] = 'C'
rs_c3, order_c3 = calc_rs(groups_c3, wedge_c3)

cip_c2 = [groups_c2[i] for i in order_c2]
cip_c3 = [groups_c3[i] for i in order_c3]
print(f"  C2 CIP优先级：{cip_c2[0]} > {cip_c2[1]} > {cip_c2[2]} > {cip_c2[3]}")
print(f"  C2 工具计算：{rs_c2} 构型")
print(f"  C3 CIP优先级：{cip_c3[0]} > {cip_c3[1]} > {cip_c3[2]} > {cip_c3[3]}")
print(f"  C3 工具计算：{rs_c3} 构型")
print(f"  目标构型：(2R, 3R)")
print(f"  验证：{'✓ C2正确' if rs_c2 == 'R' else '✗ C2错误'}，{'✓ C3正确' if rs_c3 == 'R' else '✗ C3错误'}")

# ---- 额外验证：CIP递归比较（关键科学正确性）----
print("\n" + "=" * 60)
print("【额外验证】CIP递归多层比较（修复的核心科学问题）")
print("=" * 60)
test_pairs = [
    ('C2H5', 'CH3', '乙基 > 甲基'),
    ('nPr', 'C2H5', '正丙基 > 乙基'),
    ('iPr', 'C2H5', '异丙基 > 乙基'),
    ('Ph', 'iPr', '苯基 > 异丙基'),
    ('Ph', 'Vinyl', '苯基 > 乙烯基'),
    ('iPr', 'Vinyl', '异丙基 > 乙烯基'),
    ('COOH', 'CHO', '羧基 > 醛基'),
    ('CHO', 'CH2OH', '醛基 > 羟甲基'),
    ('CH2OH', 'C2H5', '羟甲基 > 乙基'),
    ('OCH3', 'OH', '甲氧基 vs 羟基（锚点O相同，比较下一层）'),
    ('OH', 'NH2', '羟基 > 氨基'),
    ('Cl', 'OH', '氯 > 羟基'),
]
for a, b, desc in test_pairs:
    cmp_result = compare_cip(a, b)
    result = f"{a} > {b}" if cmp_result > 0 else (f"{a} < {b}" if cmp_result < 0 else f"{a} = {b}")
    print(f"  {desc}: {result} {'✓' if cmp_result > 0 else '✗'}")

print("\n" + "=" * 60)
print("验证完成")
print("=" * 60)
