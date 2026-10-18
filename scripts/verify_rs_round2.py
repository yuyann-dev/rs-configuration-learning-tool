# -*- coding: utf-8 -*-
"""
315化学(农) RS构型真题验证 - 第二轮
选取三道不同的典型题目，独立复现HTML中的calcRS和compareCIP逻辑
"""
import math

# ============ CIP规则（与HTML一致）============
GROUP_ATOM = {
    'H': 'H', 'D': 'H', 'CH3': 'C', 'C2H5': 'C', 'nPr': 'C', 'iPr': 'C',
    'nBu': 'C', 'tBu': 'C', 'Vinyl': 'C', 'Ethynyl': 'C', 'Ph': 'C',
    'benzyl': 'C', 'OH': 'O', 'OCH3': 'O', 'OEt': 'O', 'OPh': 'O',
    'NH2': 'N', 'NHCH3': 'N', 'N(CH3)2': 'N', 'NO2': 'N',
    'F': 'F', 'Cl': 'Cl', 'Br': 'Br', 'I': 'I',
    'CHO': 'C', 'COCH3': 'C', 'COOH': 'C', 'COOCH3': 'C',
    'CONH2': 'C', 'CN': 'C', 'CH2OH': 'C', 'CH2NH2': 'C',
    'CH2Cl': 'C', 'CH2Br': 'C', 'SH': 'S', 'SCH3': 'S',
    'SO3H': 'S',
}

ATOMIC_NUM = {'H':1,'C':6,'N':7,'O':8,'F':9,'P':15,'S':16,'Cl':17,'Br':35,'I':53}

# 双键/三键按CIP规则重复幻影原子
CIP_CHILDREN = {
    'H': [], 'D': [], 'F': [], 'Cl': [], 'Br': [], 'I': [],
    'CH3': ['H','H','H'],
    'C2H5': ['CH3','H','H'],
    'nPr': ['C2H5','H','H'],
    'iPr': ['CH3','CH3','H'],
    'nBu': ['nPr','H','H'],
    'tBu': ['CH3','CH3','CH3'],
    'Vinyl': ['C','C','H'],
    'Ethynyl': ['C','C','C'],
    'Ph': ['PhC','PhC','H'],
    'PhC': ['PhC','PhCH','H'],
    'PhCH': ['PhC','H'],
    'benzyl': ['Ph','H','H'],
    'OH': ['H'],
    'OCH3': ['CH3'],
    'OEt': ['C2H5'],
    'OPh': ['Ph'],
    'NH2': ['H','H'],
    'NHCH3': ['CH3','H'],
    'N(CH3)2': ['CH3','CH3'],
    'NO2': ['O','O','O'],
    'CHO': ['O','O','H'],
    'COCH3': ['O','O','CH3'],
    'COOH': ['O','O','OH'],
    'COOCH3': ['O','O','OCH3'],
    'CONH2': ['O','O','NH2'],
    'CN': ['N','N','N'],
    'CH2OH': ['OH','H','H'],
    'CH2NH2': ['NH2','H','H'],
    'CH2Cl': ['Cl','H','H'],
    'CH2Br': ['Br','H','H'],
    'SH': ['H'],
    'SCH3': ['CH3'],
    'SO3H': ['O','O','OH'],
}

def cip_atom_num(key):
    if key in GROUP_ATOM:
        return ATOMIC_NUM.get(GROUP_ATOM[key], 6)
    if key == 'PhC' or key == 'PhCH':
        return 6
    return ATOMIC_NUM.get(key, 6)

def build_cip_layers(group, max_depth=8):
    layers = []
    current = [group]
    visited = set()
    for _ in range(max_depth):
        layer = sorted([cip_atom_num(g) for g in current], reverse=True)
        layers.append(layer)
        next_level = []
        for g in current:
            if g in visited:
                continue
            visited.add(g)
            children = CIP_CHILDREN.get(g, [])
            next_level.extend(children)
        if not next_level:
            break
        current = next_level
    return layers

def compare_cip(a, b):
    an = cip_atom_num(a)
    bn = cip_atom_num(b)
    if an != bn:
        return an - bn
    la = build_cip_layers(a)
    lb = build_cip_layers(b)
    max_l = max(len(la), len(lb))
    for i in range(max_l):
        ra = la[i] if i < len(la) else []
        rb = lb[i] if i < len(lb) else []
        max_len = max(len(ra), len(rb))
        for j in range(max_len):
            va = ra[j] if j < len(ra) else 0
            vb = rb[j] if j < len(rb) else 0
            if va != vb:
                return va - vb
    return 0

# ============ 四面体坐标（与HTML一致）============
def get_tetra_coords(groups, wedge):
    base = [[0,1.25,-0.45],[1.15,0,0.35],[0,-1.25,-0.45],[-1.15,0,0.35]]
    coords = []
    for i in range(4):
        x, y, z = base[i]
        if wedge[i] == 'up': z += 0.85
        elif wedge[i] == 'down': z -= 0.85
        coords.append([x, y, z])
    return coords

def normalize(v):
    l = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    return [v[0]/l, v[1]/l, v[2]/l] if l > 0 else [0,0,0]

def rotate_point(v, axis, angle):
    x, y, z = v
    ux, uy, uz = axis
    cos_a = math.cos(angle)
    sin_a = math.sin(angle)
    dot = ux*x + uy*y + uz*z
    rx = x*cos_a + (uy*z - uz*y)*sin_a + ux*dot*(1-cos_a)
    ry = y*cos_a + (uz*x - ux*z)*sin_a + uy*dot*(1-cos_a)
    rz = z*cos_a + (ux*y - uy*x)*sin_a + uz*dot*(1-cos_a)
    return [rx, ry, rz]

def calc_rs(groups, wedge):
    coords = get_tetra_coords(groups, wedge)
    idx = sorted(range(4), key=lambda i: sum(1 for j in range(4) if compare_cip(groups[i], groups[j]) > 0), reverse=True)
    ordered = [coords[i] for i in idx]
    lowest = ordered[3]
    n = normalize(lowest)
    target = [0, 0, -1]
    axis = [n[1]*target[2]-n[2]*target[1], n[2]*target[0]-n[0]*target[2], n[0]*target[1]-n[1]*target[0]]
    al = math.sqrt(axis[0]**2 + axis[1]**2 + axis[2]**2)
    if al < 1e-6:
        axis = [1, 0, 0]
    else:
        axis = [axis[0]/al, axis[1]/al, axis[2]/al]
    angle = math.acos(max(-1, min(1, n[0]*target[0] + n[1]*target[1] + n[2]*target[2])))
    r1 = rotate_point(ordered[0], axis, angle)
    r2 = rotate_point(ordered[1], axis, angle)
    r3 = rotate_point(ordered[2], axis, angle)
    cross2d = (r2[0]-r1[0])*(r3[1]-r1[1]) - (r2[1]-r1[1])*(r3[0]-r1[0])
    return ('S' if cross2d > 0 else 'R'), idx

# ============ 验证题目 ============
print("=" * 60)
print("315化学(农) RS构型真题验证 - 第二轮（三道新题）")
print("=" * 60)

# ---- 基础题：(S)-2-氯丁烷 ----
print("\n【基础题】(S)-2-氯丁烷")
print("结构：CH3-CH(Cl)-CH2-CH3")
print("手性碳C2连接：Cl, C2H5, CH3, H")
groups1 = ['C2H5', 'H', 'CH3', 'Cl']
wedge1 = ['plain', 'down', 'plain', 'up']
rs1, order1 = calc_rs(groups1, wedge1)
cip1 = [groups1[i] for i in order1]
print(f"  CIP优先级：{cip1[0]} > {cip1[1]} > {cip1[2]} > {cip1[3]}")
print(f"  工具计算结果：{rs1} 构型")
print(f"  正确答案：S 构型")
print(f"  验证：{'✓ 正确' if rs1 == 'S' else '✗ 错误'}")

# ---- 中档题：D-甘油醛 Fischer投影式 ----
print("\n【中档题】D-甘油醛 Fischer投影式判断")
print("Fischer投影式：CHO在上, OH在右, CH2OH在下, H在左")
print("（横前竖后，H在横键朝前）")
groups2 = ['CHO', 'OH', 'CH2OH', 'H']
wedge2 = ['down', 'up', 'down', 'up']
rs2, order2 = calc_rs(groups2, wedge2)
cip2 = [groups2[i] for i in order2]
print(f"  CIP优先级：{cip2[0]} > {cip2[1]} > {cip2[2]} > {cip2[3]}")
print(f"  工具计算结果：{rs2} 构型")
print(f"  正确答案：R 构型（D-甘油醛）")
print(f"  验证：{'✓ 正确' if rs2 == 'R' else '✗ 错误'}")

# ---- 难题：L-苏氨酸 (2S,3R) ----
print("\n【难题】L-苏氨酸 (2S,3R)-2-氨基-3-羟基丁酸")
print("结构：CH3-CH(OH)-CH(NH2)-COOH")
# C2 Fischer: COOH上, NH2左(前), C3下, H右(前)
# 转楔形键[top,right,bottom,left]: COOH, H, C3, NH2; wedge: down, up, down, up
groups_c2 = ['COOH', 'H', '_frag_c3', 'NH2']
wedge_c2 = ['down', 'up', 'down', 'up']
CIP_CHILDREN['_frag_c3'] = ['OH', 'CH3', 'H']
GROUP_ATOM['_frag_c3'] = 'C'
rs_c2, order_c2 = calc_rs(groups_c2, wedge_c2)

# C3 Fischer: C2上, OH右(前), CH3下, H左(前)
groups_c3 = ['_frag_c2', 'OH', 'CH3', 'H']
wedge_c3 = ['down', 'up', 'down', 'up']
CIP_CHILDREN['_frag_c2'] = ['NH2', 'COOH', 'H']
GROUP_ATOM['_frag_c2'] = 'C'
rs_c3, order_c3 = calc_rs(groups_c3, wedge_c3)

cip_c2 = [groups_c2[i] for i in order_c2]
cip_c3 = [groups_c3[i] for i in order_c3]
print(f"  C2 CIP优先级：{cip_c2[0]} > {cip_c2[1]} > {cip_c2[2]} > {cip_c2[3]}")
print(f"  C2 工具计算：{rs_c2} 构型")
print(f"  C3 CIP优先级：{cip_c3[0]} > {cip_c3[1]} > {cip_c3[2]} > {cip_c3[3]}")
print(f"  C3 工具计算：{rs_c3} 构型")
print(f"  目标构型：(2S, 3R)")
print(f"  验证：{'✓ C2正确' if rs_c2 == 'S' else '✗ C2错误'}，{'✓ C3正确' if rs_c3 == 'R' else '✗ C3错误'}")

# ============ 额外CIP验证 ============
print("\n" + "=" * 60)
print("【额外验证】CIP递归多层比较")
print("=" * 60)
tests = [
    ('Br', 'Cl', '溴 > 氯'),
    ('Cl', 'OH', '氯 > 羟基'),
    ('OH', 'NH2', '羟基 > 氨基'),
    ('NH2', 'COOH', '氨基 > 羧基（锚点N>C）'),
    ('COOH', 'CHO', '羧基 > 醛基'),
    ('CHO', 'CH2OH', '醛基 > 羟甲基'),
    ('CH2OH', 'C2H5', '羟甲基 > 乙基'),
    ('iPr', 'C2H5', '异丙基 > 乙基'),
    ('tBu', 'iPr', '叔丁基 > 异丙基'),
    ('Ph', 'Vinyl', '苯基 > 乙烯基'),
    ('Ethynyl', 'Vinyl', '乙炔基 > 乙烯基'),
    ('CN', 'COOH', '氰基 vs 羧基'),
]
for a, b, desc in tests:
    r = compare_cip(a, b)
    ok = '✓' if r > 0 else '✗'
    print(f"  {desc}: {a} {'>' if r>0 else '<'} {b} {ok}")

print("\n" + "=" * 60)
print("验证完成")
print("=" * 60)
