# -*- coding: utf-8 -*-
"""验证3道真题：预设库 + 自定义构建器双路径"""
import math

ATOMIC_NUM = {'H':1,'D':2,'C':6,'N':7,'O':8,'F':9,'P':15,'S':16,'Cl':17,'Br':35,'I':53}
def atomic_num(e): return ATOMIC_NUM.get(e, 0)

CIP_CHILDREN = {
    'H':[],'D':[],'F':[],'Cl':[],'Br':[],'I':[],
    'CH3':['H','H','H'],'C2H5':['CH3','H','H'],'nPr':['C2H5','H','H'],
    'iPr':['CH3','CH3','H'],'tBu':['CH3','CH3','CH3'],
    'Vinyl':['C','C','H'],'Ethynyl':['C','C','C'],'Allyl':['Vinyl','H','H'],
    'Ph':['PhC','PhC','H'],'PhC':['PhC','PhCH','H'],'PhCH':['PhC','H'],'Bn':['Ph','H','H'],
    'OH':['H'],'OCH3':['CH3'],'OAc':['COCH3'],
    'CHO':['O','O','H'],'COCH3':['O','O','CH3'],'COOH':['O','O','OH'],
    'CH2OH':['OH','H','H'],'CH2NH2':['NH2','H','H'],'CN':['N','N','N'],
    'NH2':['H','H'],'NO2':['O','O','O'],
    'SH':['H'],'SCH3':['CH3'],
    'CHOHCOOH':['OH','COOH','H'],'CHOHCH3':['OH','CH3','H'],
    'CHClCH3':['Cl','CH3','H'],'CHNH2COOH':['NH2','COOH','H'],
}
GROUP_ATOM = {
    'H':'H','D':'H','F':'F','Cl':'Cl','Br':'Br','I':'I',
    'CH3':'C','C2H5':'C','nPr':'C','iPr':'C','tBu':'C',
    'Vinyl':'C','Ethynyl':'C','Allyl':'C','Ph':'C','Bn':'C',
    'OH':'O','OCH3':'O','OAc':'O','CHO':'C','COCH3':'C','COOH':'C',
    'CH2OH':'C','CH2NH2':'C','CN':'C','NH2':'N','NO2':'N','SH':'S','SCH3':'S',
    'CHOHCOOH':'C','CHOHCH3':'C','CHClCH3':'C','CHNH2COOH':'C',
}
def cip_atom_num(k):
    if k in GROUP_ATOM: return atomic_num(GROUP_ATOM[k])
    if k.startswith('=') or k.startswith('#'): return atomic_num(k[1])
    if k == 'PhC': return 6
    return atomic_num(k) or 6

cip_cache = {}
def build_cip_layers(key, depth=0, visited=None):
    if depth > 8: return []
    ck = key + '_' + str(depth)
    if ck in cip_cache: return cip_cache[ck]
    if not visited: visited = set()
    if key in visited: return []
    visited.add(key)
    children = CIP_CHILDREN.get(key, [])
    if not children: cip_cache[ck] = []; return []
    layer = sorted([cip_atom_num(c) for c in children], reverse=True)
    next_layers = []
    for c in children:
        sub = build_cip_layers(c, depth+1, set(visited))
        for i, s in enumerate(sub):
            if i >= len(next_layers): next_layers.append([])
            next_layers[i].extend(s)
    result = [layer] + [sorted(l, reverse=True) for l in next_layers]
    cip_cache[ck] = result
    return result

def compare_cip(a, b):
    an, bn = cip_atom_num(a), cip_atom_num(b)
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

def normalize(v):
    l = math.hypot(v[0], v[1], v[2]) or 1
    return [v[0]/l, v[1]/l, v[2]/l]

def rotate_point(v, axis, angle):
    c, s = math.cos(angle), math.sin(angle)
    x, y, z = v; ux, uy, uz = axis
    return [
        x*(c+ux*ux*(1-c)) + y*(ux*uy*(1-c)-uz*s) + z*(ux*uz*(1-c)+uy*s),
        x*(uy*ux*(1-c)+uz*s) + y*(c+uy*uy*(1-c)) + z*(uy*uz*(1-c)-ux*s),
        x*(uz*ux*(1-c)-uy*s) + y*(uz*uy*(1-c)+ux*s) + z*(c+uz*uz*(1-c))
    ]

def calc_rs(center):
    base = [[0,1.25,-0.45],[1.15,0,0.35],[0,-1.25,-0.45],[-1.15,0,0.35]]
    coords = []
    for i, w in enumerate(center['wedge']):
        x, y, z = base[i]
        if w == 'up': z += 0.85
        elif w == 'down': z -= 0.85
        coords.append([x, y, z])
    idx = [0,1,2,3]
    for i in range(4):
        for j in range(i+1, 4):
            if compare_cip(center['groups'][idx[j]], center['groups'][idx[i]]) > 0:
                idx[i], idx[j] = idx[j], idx[i]
    center['cipOrder'] = idx
    ordered = [coords[i] for i in idx]
    v4 = ordered[3]; n = normalize(v4); target = [0,0,-1]
    axis = [n[1]*target[2]-n[2]*target[1], n[2]*target[0]-n[0]*target[2], n[0]*target[1]-n[1]*target[0]]
    al = math.hypot(*axis)
    angle = math.acos(max(-1, min(1, n[0]*target[0]+n[1]*target[1]+n[2]*target[2])))
    if al < 1e-6:
        r1, r2, r3 = ordered[0], ordered[1], ordered[2]
    else:
        axis = [a/al for a in axis]
        r1 = rotate_point(ordered[0], axis, angle)
        r2 = rotate_point(ordered[1], axis, angle)
        r3 = rotate_point(ordered[2], axis, angle)
    cross2d = (r2[0]-r1[0])*(r3[1]-r1[1]) - (r2[1]-r1[1])*(r3[0]-r1[0])
    return 'S' if cross2d > 0 else 'R'

passed = 0; total = 0
def check(name, got, expect):
    global passed, total
    total += 1
    ok = got == expect
    if ok: passed += 1
    print(f"  {name}: {got} (期望 {expect}) {'✓' if ok else '✗ FAIL'}")

print("=" * 55)
print("315化学(农)真题验证 — 预设库 + 自定义构建器")
print("=" * 55)

# 题目1：2-氯丁烷（基础卤代烃）
print("\n【题1】2-氯丁烷 CH3-CH(Cl)-CH2-CH3")
print("  CIP: Cl(17) > C2H5 > CH3 > H")
s_preset = {'groups':['C2H5','H','CH3','Cl'], 'wedge':['plain','down','plain','up']}
check("预设库 S-2-氯丁烷", calc_rs(s_preset), 'S')
r_builder = {'groups':['C2H5','Cl','CH3','H'], 'wedge':['plain','up','plain','down']}
check("自定义构建 R-2-氯丁烷", calc_rs(r_builder), 'R')

# 题目2：甘油醛（中档，糖类）
print("\n【题2】甘油醛 CHO-CH(OH)-CH2OH")
print("  CIP: OH > CHO > CH2OH > H")
r_preset = {'groups':['CHO','OH','CH2OH','H'], 'wedge':['plain','up','plain','down']}
check("预设库 R-甘油醛(D型)", calc_rs(r_preset), 'R')
s_builder = {'groups':['CHO','H','CH2OH','OH'], 'wedge':['plain','down','plain','up']}
check("自定义构建 S-甘油醛(L型)", calc_rs(s_builder), 'S')

# 题目3：酒石酸（难题，双手性中心）
print("\n【题3】酒石酸 HOOC-CH(OH)-CH(OH)-COOH")
print("  CIP: OH > COOH > CH(OH)COOH > H")
rr_c1 = {'groups':['COOH','OH','CHOHCOOH','H'], 'wedge':['plain','up','plain','down']}
rr_c2 = {'groups':['COOH','OH','CHOHCOOH','H'], 'wedge':['plain','up','plain','down']}
check("预设库 (R,R)-酒石酸 C2", calc_rs(rr_c1), 'R')
check("预设库 (R,R)-酒石酸 C3", calc_rs(rr_c2), 'R')
# 自定义构建器双手性中心
b_c1 = {'groups':['COOH','OH','CHOHCOOH','H'], 'wedge':['plain','up','plain','down']}
b_c2 = {'groups':['COOH','OH','CHOHCOOH','H'], 'wedge':['plain','up','plain','down']}
check("自定义构建 (R,R) C2", calc_rs(b_c1), 'R')
check("自定义构建 (R,R) C3", calc_rs(b_c2), 'R')

# 附加：内消旋酒石酸
print("\n【附加】meso-酒石酸 (2R,3S)")
meso_c1 = {'groups':['COOH','OH','CHOHCOOH','H'], 'wedge':['plain','up','plain','down']}
meso_c2 = {'groups':['COOH','H','CHOHCOOH','OH'], 'wedge':['plain','down','plain','up']}
check("meso C2", calc_rs(meso_c1), 'R')
check("meso C3", calc_rs(meso_c2), 'S')

print(f"\n{'='*55}")
print(f"结果: {passed}/{total} 通过")
print("=" * 55)
