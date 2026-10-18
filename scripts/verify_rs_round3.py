# -*- coding: utf-8 -*-
"""
315化学(农) RS构型真题验证 - 第三轮
全部使用自定义结构模式（模拟二维构建器输入），不使用预设化合物
三道题：基础(含异丙基+卤素)、中档(内消旋酒石酸)、难题(异亮氨酸双手性)
"""
import math

# ============ CIP规则（与HTML完全一致）============
GROUP_ATOM = {
    'H': 'H', 'D': 'H', 'CH3': 'C', 'C2H5': 'C', 'nPr': 'C', 'iPr': 'C',
    'nBu': 'C', 'tBu': 'C', 'Vinyl': 'C', 'Ethynyl': 'C', 'Ph': 'C',
    'benzyl': 'C', 'OH': 'O', 'OCH3': 'O', 'OEt': 'O', 'OPh': 'O',
    'NH2': 'N', 'NHCH3': 'N', 'N(CH3)2': 'N', 'NO2': 'N',
    'F': 'F', 'Cl': 'Cl', 'Br': 'Br', 'I': 'I',
    'CHO': 'C', 'COCH3': 'C', 'COOH': 'C', 'COOCH3': 'C',
    'CONH2': 'C', 'CN': 'C', 'CH2OH': 'C', 'CH2NH2': 'C',
    'CH2Cl': 'C', 'CH2Br': 'C', 'SH': 'S', 'SCH3': 'S', 'SO3H': 'S',
}

ATOMIC_NUM = {'H':1,'C':6,'N':7,'O':8,'F':9,'P':15,'S':16,'Cl':17,'Br':35,'I':53}

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
    if key in ('PhC', 'PhCH'):
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
def get_tetra_coords(wedge):
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
    coords = get_tetra_coords(wedge)
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

# ============ 自定义构建器输入验证 ============
print("=" * 65)
print("315化学(农) RS构型验证 - 第三轮（纯自定义结构模式）")
print("=" * 65)

# ===== 基础题：(R)-2-溴-3-甲基丁烷 =====
print("\n【基础题】(R)-2-溴-3-甲基丁烷")
print("自定义构建器输入（单手性模式）：")
print("  上基团=iPr, 右基团=Br(实楔↑), 下基团=CH3, 左基团=H(虚楔↓)")
groups1 = ['iPr', 'Br', 'CH3', 'H']
wedge1 = ['plain', 'up', 'plain', 'down']
rs1, order1 = calc_rs(groups1, wedge1)
cip1 = [groups1[i] for i in order1]
print(f"  CIP优先级：{cip1[0]} > {cip1[1]} > {cip1[2]} > {cip1[3]}")
print(f"  工具计算：{rs1} 构型")
print(f"  手动验证：Br>iPr>CH3>H, H在左(后), Br(右)→iPr(上)→CH3(下): 顺时针 → R")
print(f"  正确答案：R 构型")
print(f"  验证：{'✓ 正确' if rs1 == 'R' else '✗ 错误'}")

# ===== 中档题：内消旋酒石酸 (2R,3S) =====
print("\n【中档题】内消旋酒石酸 (2R,3S)-2,3-二羟基丁二酸")
print("自定义构建器输入（双手性模式）：")
print("  上端基=COOH, 下端基=COOH")
print("  C1: 右=OH(实楔↑), 左=H(虚楔↓)")
print("  C2: 右=OH(实楔↑), 左=H(虚楔↓)")
# C2(上端手性碳): top=COOH, right=OH(up), bottom=C3片段, left=H(down)
groups_c2 = ['COOH', 'OH', '_frag_meso_c3', 'H']
wedge_c2 = ['plain', 'up', 'plain', 'down']
CIP_CHILDREN['_frag_meso_c3'] = ['OH', 'COOH', 'H']
GROUP_ATOM['_frag_meso_c3'] = 'C'
rs_c2, order_c2 = calc_rs(groups_c2, wedge_c2)

# C3(下端手性碳): top=C2片段, right=OH(up), bottom=COOH, left=H(down)
# 注意：链状双手性中心保持一致的左右朝向，构型差异来自上下基团CIP排序不同
groups_c3 = ['_frag_meso_c2', 'OH', 'COOH', 'H']
wedge_c3 = ['plain', 'up', 'plain', 'down']
CIP_CHILDREN['_frag_meso_c2'] = ['OH', 'COOH', 'H']
GROUP_ATOM['_frag_meso_c2'] = 'C'
rs_c3, order_c3 = calc_rs(groups_c3, wedge_c3)

cip_c2 = [groups_c2[i] for i in order_c2]
cip_c3 = [groups_c3[i] for i in order_c3]
print(f"  C2 CIP：{cip_c2[0]} > {cip_c2[1]} > {cip_c2[2]} > {cip_c2[3]}")
print(f"  C2 工具计算：{rs_c2} 构型")
print(f"  C3 CIP：{cip_c3[0]} > {cip_c3[1]} > {cip_c3[2]} > {cip_c3[3]}")
print(f"  C3 工具计算：{rs_c3} 构型")
print(f"  目标：(2R, 3S) 内消旋（分子有对称面，无旋光性）")
print(f"  验证：{'✓ C2正确' if rs_c2 == 'R' else '✗ C2错误'}，{'✓ C3正确' if rs_c3 == 'S' else '✗ C3错误'}")
print(f"  内消旋判断：C2和C3构型相反({'✓ 是内消旋' if rs_c2 != rs_c3 else '✗ 不是内消旋'})")

# ===== 难题：L-异亮氨酸 (2S,3S) =====
print("\n【难题】L-异亮氨酸 (2S,3S)-2-氨基-3-甲基戊酸")
print("自定义构建器输入（双手性模式）：")
print("  上端基=COOH, 下端基=C2H5")
print("  C1(α碳): 右=H(虚楔↓), 左=NH2(实楔↑)")
print("  C2(β碳): 右=H(虚楔↓), 左=CH3(实楔↑)")
# C2(α碳): top=COOH, right=H(down), bottom=C3片段, left=NH2(up)
groups_ile_c2 = ['COOH', 'H', '_frag_ile_c3', 'NH2']
wedge_ile_c2 = ['plain', 'down', 'plain', 'up']
CIP_CHILDREN['_frag_ile_c3'] = ['CH3', 'C2H5', 'H']
GROUP_ATOM['_frag_ile_c3'] = 'C'
rs_ile_c2, order_ile_c2 = calc_rs(groups_ile_c2, wedge_ile_c2)

# C3(β碳): top=C2片段, right=H(down), bottom=C2H5, left=CH3(up)
groups_ile_c3 = ['_frag_ile_c2', 'H', 'C2H5', 'CH3']
wedge_ile_c3 = ['plain', 'down', 'plain', 'up']
CIP_CHILDREN['_frag_ile_c2'] = ['NH2', 'COOH', 'H']
GROUP_ATOM['_frag_ile_c2'] = 'C'
rs_ile_c3, order_ile_c3 = calc_rs(groups_ile_c3, wedge_ile_c3)

cip_ile_c2 = [groups_ile_c2[i] for i in order_ile_c2]
cip_ile_c3 = [groups_ile_c3[i] for i in order_ile_c3]
print(f"  C2(α) CIP：{cip_ile_c2[0]} > {cip_ile_c2[1]} > {cip_ile_c2[2]} > {cip_ile_c2[3]}")
print(f"  C2(α) 工具计算：{rs_ile_c2} 构型")
print(f"  C3(β) CIP：{cip_ile_c3[0]} > {cip_ile_c3[1]} > {cip_ile_c3[2]} > {cip_ile_c3[3]}")
print(f"  C3(β) 工具计算：{rs_ile_c3} 构型")
print(f"  目标：(2S, 3S) L-异亮氨酸")
print(f"  验证：{'✓ C2正确' if rs_ile_c2 == 'S' else '✗ C2错误'}，{'✓ C3正确' if rs_ile_c3 == 'S' else '✗ C3错误'}")

# ============ 额外CIP边界验证 ============
print("\n" + "=" * 65)
print("【额外验证】CIP边界情况比较")
print("=" * 65)
tests = [
    ('Br', 'iPr', '溴 > 异丙基（锚点Br>C）'),
    ('iPr', 'C2H5', '异丙基 > 乙基（第二层C>C vs C>H）'),
    ('iPr', 'nPr', '异丙基 > 正丙基（第二层C,C vs C,H）'),
    ('COOH', 'CN', '羧基 > 氰基（O>N）'),
    ('CHO', 'Vinyl', '醛基 > 乙烯基（O>C）'),
    ('Ethynyl', 'iPr', '乙炔基 vs 异丙基'),
    ('Ph', 'tBu', '苯基 vs 叔丁基'),
    ('NO2', 'COOH', '硝基 > 羧基（锚点N>C？不，N<C...）'),
]
for a, b, desc in tests:
    r = compare_cip(a, b)
    symbol = '>' if r > 0 else ('<' if r < 0 else '=')
    print(f"  {desc}: {a} {symbol} {b}")

print("\n" + "=" * 65)
print("验证完成")
print("=" * 65)
