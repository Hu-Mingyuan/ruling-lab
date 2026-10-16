"""Exact independent geometric check of the connected S7 nonprimitive eye.

Python 3.10+, no third-party imports and no ruling/tropical engine imports.
Run: python verify_nonprimitive_witness.py --output witness-check.json
The small JSON beside this file contains four universal-cover polygonal paths
and the original four straight geodesics, not an answer table or engine dump.
Coordinates are scaled by 7063. The checker verifies complete single coverage
of every geodesic, eye geometry, all six switches, type-B companion orders,
all basepoint translates, connectedness, and the genuine nonprimitive period.
Nonswitch sector geometry is recorded explicitly. Both opposite-codirection
models are treated as overlaps, as in the manuscript's crossing-types remark.
"""
import argparse
import json
from collections import Counter, defaultdict
from fractions import Fraction as F
from math import gcd
from pathlib import Path


def sign(x):
    return (x > 0) - (x < 0)


def bezout(a,b):
    if b==0:return (1 if a>0 else -1),0
    q,r=divmod(a,b); x,y=bezout(b,r)
    return y,x-q*y


def verify(data):
    L=data['period']; bx,by=map(F,data['base_point'])
    paths=data['paths']; assert [q['id'] for q in paths]==list(range(4))
    disks=data['disk_paths']; lower,upper=data['annular_pair']
    assert data['annular_relative_deck']==[0,0]
    geodesics=data['geodesics']; intervals=defaultdict(list); nodes=defaultdict(list)
    for path in paths:
        pts=path['points']; hx,hy=path['homology']; path['segments']=[]
        path['visits']=[]
        for k,(x,y) in enumerate(pts):
            nx,ny=pts[k+1] if k+1<len(pts) else (pts[0][0]+hx*L,pts[0][1]+hy*L)
            px,py=pts[k-1] if k else (pts[-1][0]-hx*L,pts[-1][1]-hy*L)
            dx,dy=nx-x,ny-y; assert dx
            matches=[]
            for geo in geodesics:
                gx,gy=geo['direction']; qx,qy=geo['point']
                if dx*gy-dy*gx==0 and (-gy*(x-qx)+gx*(y-qy))%L==0:
                    matches.append(geo)
            assert len(matches)==1
            geo=matches[0]; gx,gy=geo['direction']; qx,qy=geo['point']
            assert gcd(gx,gy)==1 and gx>0
            assert geo['coorientation_y_sign']*sign(dx)==path['sign']
            aa,bb=bezout(gx,gy); assert aa*gx+bb*gy==1
            sx,sy=(x,y) if dx>0 else (nx,ny)
            start=F(aa*(sx-qx)+bb*(sy-qy),L)%1
            length=F(abs(dx),L*gx); assert 0<length<1
            end=start+length
            intervals[geo['id']].append((start,min(end,F(1))))
            if end>1:intervals[geo['id']].append((F(0),end-1))
            path['segments'].append((x,y,nx,ny))
            visit=dict(path=path['id'],point=(x,y),incoming=(x-px,y-py),outgoing=(dx,dy))
            nodes[(x%L,y%L)].append(visit); path['visits'].append(visit)
        path['area2']=sum(x*v-y*u for x,y,u,v in path['segments'])
    # These intervals prove that the four paths use every open geodesic arc
    # exactly once; no sampling or numerical proximity tolerance is used.
    for geo in geodesics:
        cur=F(0)
        for start,end in sorted(intervals[geo['id']]):
            assert start==cur and end>start;cur=end
        assert cur==1
    assert len(nodes)==28 and all(len(q)==2 for q in nodes.values())

    def hits(i,x):
        path=paths[i]; hx,hy=path['homology']; out=set()
        if hx<0:hx,hy=-hx,-hy
        for x0,y0,x1,y1 in path['segments']:
            if x0>x1:x0,y0,x1,y1=x1,y1,x0,y0
            assert x0<x1
            ns=[0] if hx==0 else range((x-x1)//(hx*L),(x-x1)//(hx*L)+3)
            for n in ns:
                u=x0+n*hx*L;v=x1+n*hx*L
                if u<=x<=v:out.add(F(y0+n*hy*L)+(x-u)*F(y1-y0,x1-x0))
        return sorted(out)

    for i in disks:
        path=paths[i]; assert path['homology']==[0,0]
        assert path['color']==path['sign']*sign(path['area2'])
        xs=sorted({x for x,y in path['points']})
        assert len(hits(i,F(xs[0])))==len(hits(i,F(xs[-1])))==1
        assert all(len(hits(i,F(x)))==2 for x in xs[1:-1])
        assert all(len(hits(i,F(x+y,2)))==2 for x,y in zip(xs,xs[1:]))
    h=paths[lower]['homology'];assert h==paths[upper]['homology']==[2,-2]
    assert paths[lower]['sign']==-1 and paths[upper]['sign']==1
    assert paths[lower]['color']==paths[upper]['color']==-1
    assert all(x1>x0 for i in (lower,upper) for x0,y0,x1,y1 in paths[i]['segments'])
    breaks=sorted({F(0),F(2*L)}|{F(x%(2*L)) for i in (lower,upper) for x,y in paths[i]['points']})
    gaps=[]
    for x in breaks:
        low=hits(lower,x);high=hits(upper,x);assert len(low)==len(high)==1
        gaps.append(high[0]-low[0])
    assert min(gaps)>0  # Difference is affine between common breakpoints.
    half_period={i:hits(i,F(L))[0]+L-hits(i,F(0))[0] for i in (lower,upper)}
    assert all(value!=0 for value in half_period.values())

    disk_multiplicities=[]
    for i in disks:
        xs=[x for x,y in paths[i]['points']];ys=[y for x,y in paths[i]['points']];count=0
        for k in range((bx-max(xs))//L-1,(bx-min(xs))//L+2):
            hh=hits(i,bx-k*L)
            if len(hh)!=2:continue
            for ell in range((by-max(ys))//L-1,(by-min(ys))//L+2):
                count+=hh[0]<by-ell*L<hh[1]
        disk_multiplicities.append(count)
    annular_multiplicities=[]
    for residue in range(2):
        low=hits(lower,bx-residue*L)[0];high=hits(upper,bx-residue*L)[0]
        annular_multiplicities.append(sum(low<by-ell*L<high
            for ell in range((by-high)//L-1,(by-low)//L+2)))
    balance=sum(paths[i]['color']*n for i,n in zip(disks,disk_multiplicities))-sum(annular_multiplicities)
    assert balance==0

    eye={disks[0]:0,disks[1]:1,lower:2,upper:2}
    interior={i:sign(paths[i]['area2']) for i in disks};interior.update({lower:1,upper:-1})
    local=[];edges=[];counts=Counter();eps=F(1,10)
    for vertex,visits in sorted(nodes.items()):
        x,y=vertex;info=[]
        for visit in visits:
            i=visit['path'];ix,iy=visit['incoming'];ox,oy=visit['outgoing']
            cross=ix*oy-iy*ox
            row=dict(path=i,color=paths[i]['color'],turn=sign(cross),convex=cross*interior[i]>0,
                     tip=ix*ox<0,near=[],far=[])
            shiftx=visit['point'][0]-x;shifty=visit['point'][1]-y
            for e in (-eps,eps):
                rays=[(-ix,-iy),(ox,oy)]
                visible=[v for v in rays if v[0]*e>0]
                if row['tip']:continue
                assert len(visible)==1
                vx,vy=visible[0];near=F(y)+e*F(vy,vx)
                if i in disks:
                    companions=[v-shifty for v in hits(i,x+e+shiftx) if v-shifty!=near]
                else:
                    companion=upper if i==lower else lower
                    companions=[v-shifty for v in hits(companion,x+e+shiftx)]
                assert len(companions)==1
                row['near'].append(near);row['far'].append(companions[0])
            info.append(row)
        switched=bool(info[0]['turn']);assert switched==bool(info[1]['turn'])
        aa,bb=info
        sectors=None
        if switched:
            assert eye[aa['path']]!=eye[bb['path']]
            if aa['tip'] or bb['tip']:
                assert aa['tip'] and bb['tip'] and aa['color']!=bb['color'] and aa['convex'] and bb['convex']
                kind='horizontal A'
            elif aa['color']!=bb['color']:
                assert aa['convex'] and bb['convex'];kind='vertical A'
            else:
                assert aa['convex']!=bb['convex'];kind='B'
                assert all((aa['near'][k]-bb['near'][k])*(aa['far'][k]-bb['far'][k])>0 for k in (0,1))
            edges.append((eye[aa['path']],eye[bb['path']]))
        else:
            same_color=aa['color']==bb['color']
            sidea=sign(aa['far'][0]-aa['near'][0]);sideb=sign(bb['far'][0]-bb['near'][0])
            assert sidea and sideb
            kind='overlap' if same_color!=(sidea==sideb) else 'departure/return'
            # A nonswitch has two unchanged straight strands, and each incident
            # eye lies strictly on the prescribed side on both nearby slices.
            assert all(sign(row['far'][k]-row['near'][k])==sign(row['far'][0]-row['near'][0])
                       for row in info for k in (0,1))
            slopes=[F(v['outgoing'][1],v['outgoing'][0]) for v in visits]
            middle=sum(slopes)/2
            samples={'N':(F(x),F(y)+eps),'S':(F(x),F(y)-eps),
                     'W':(F(x)-eps,F(y)-middle*eps),'E':(F(x)+eps,F(y)+middle*eps)}
            sectors={}
            for name,(xx,yy) in samples.items():
                occupied=[]
                for v,row,m in zip(visits,info,slopes):
                    i=v['path'];sx=v['point'][0]-x;sy=v['point'][1]-y
                    bounds=hits(i,xx+sx) if i in disks else [hits(lower,xx+sx)[0],hits(upper,xx+sx)[0]]
                    assert len(bounds)==2 and bounds[0]<bounds[1]
                    contained=bounds[0]-sy<yy<bounds[1]-sy
                    side=sign(row['far'][0]-row['near'][0])
                    assert contained==(side*(yy-y-m*(xx-x))>0)
                    occupied.append(contained)
                sectors[name]=occupied
        counts[kind]+=1;local.append(dict(vertex=vertex,type=kind,branches=info,sectors=sectors))
    assert counts['horizontal A']==counts['vertical A']==counts['B']==2
    connected={0}
    while True:
        more={j for i,j in edges if i in connected}|{i for i,j in edges if j in connected}
        if more<=connected:break
        connected|=more
    assert connected=={0,1,2}
    chi=len(disks)-len(edges);boundary=len(geodesics);genus=F(2-boundary-chi,2)
    assert (chi,boundary,genus)==(-4,4,1)
    return dict(verified=True,annular_period=h,annular_gap_min=min(gaps),annular_gap_max=max(gaps),
                half_period_defect=half_period,disk_base_multiplicities=disk_multiplicities,
                annular_base_multiplicities_by_residue=annular_multiplicities,balance=balance,
                local_type_counts=counts,local_checks=local,connected=True,
                euler_characteristic=chi,boundary_components=boundary,genus=int(genus),
                original_geodesic_arcs_covered_once=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('geometry',nargs='?',type=Path,
                        default=Path(__file__).with_name('nonprimitive-witness-geometry.json'))
    parser.add_argument('--output',type=Path)
    args=parser.parse_args();result=verify(json.loads(args.geometry.read_text(encoding='utf-8')))
    if args.output:args.output.write_text(json.dumps(result,indent=2,default=str)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='local_checks'},default=str))


if __name__=='__main__':main()
