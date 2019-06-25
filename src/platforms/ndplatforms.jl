export NdBSplinePlatform
"""
    NdBSplinePlatform(orders::NTuple{N,Int), types=ntuple(k=>BSplinePlatform, Val(N))) where N

Return a ProductPlatform that has B-spline platforms of type `types` and order `orders` as its elements.
"""
NdBSplinePlatform(orders::NTuple{N,Int}, types = ntuple(k->BSplinePlatform, Val(N))) where N =
    FrameFun.ProductPlatform(map((d,type)->type(d), orders, types)...)

export NdEpsBSplinePlatform
"""
    NdEpsBSplinePlatform(orders::NTuple{N,Int), types=ntuple(k=>BSplinePlatform, Val(N))) where N

Return a ProductPlatform that has B-spline platforms of type `types` and order `orders` as its elements.
"""
NdEpsBSplinePlatform(orders::NTuple{N,Int}, types = ntuple(k->EpsBSplinePlatform, Val(N))) where N =
    FrameFun.ProductPlatform(map((d,type)->type(d), orders, types)...)

export NdCDBSplinePlatform
"""
    NdCDBSplinePlatform(orders::NTuple{N,Int), types=ntuple(k=>BSplinePlatform, Val(N))) where N

Return a ProductPlatform that has B-spline platforms of type `types` and order `orders` as its elements.
"""
NdCDBSplinePlatform(orders::NTuple{N,Int}, types = ntuple(k->CDBSplinePlatform, Val(N))) where N =
    FrameFun.ProductPlatform(map((d,type)->type(d), orders, types)...)
