## Date:2015-06-23
## Author:wonjin

immutable Point
    coord::Vector{Float64}
    function Point(x, y, z)
        coord = [x, y, z]
        new(coord)
    end
end

+(p1::Point, p2::Point) = Point(p1.coord[1]+p2.coord[1], p1.coord[2]+p2.coord[2], p1.coord[3]+p2.coord[3])
-(p1::Point, p2::Point) = Point(p1.coord[1]-p2.coord[1], p1.coord[2]-p2.coord[2], p1.coord[3]-p2.coord[3])
*(p1::Point, s::Real) = Point(u*p1.coord[1], u*p1.coord[2], u*p1.coord[3])
*(s::Real, p1::Point) = Point(u*p1.coord[1], u*p1.coord[2], u*p1.coord[3])
getx(p::Point) = p.coord[1]
gety(p::Point) = p.coord[2]
getz(p::Point) = p.coord[3]

immutable Line
    p::Point
    direction::Vector{Float64}
    Line(p::Point, dir::Vector{Float64}) = new(p, dir/norm(dir))
    Line(p::Point, dir) = new(p, convert(Vector{Float64}, dir/norm(dir)))
end
function Line(p1::Point, p2::Point)
    diff = p2 - p1
    Line(p1, diff.coord)
end

function isparallel(l1::Line, l2::Line)
    if isapprox(l1.direction[1], l2.direction[1]) &&
        isapprox(l1.direction[2], l2.direction[2]) &&
        isapprox(l1.direction[3], l2.direction[3])
        return true
    end
    return false
end

function findjunction(l1::Line, l2::Line)
    if isparallel(l1, l2)
        throw(error("lines must not be parallel"))
    end
    x = (gety(l2.p) - gety(l1.p) - l2.direction[2]/l2.direction[1]*getx(l2.p)
    + l1.direction[2]/l1.direction[1]*getx(l1.p)) / (l1.direction[2]/l1.direction[1]
    - l2.direction[2]/l2.direction[1])
    y = l1.direction[2]/l1.direction[1]*(x-getx(l1.p)) + gety(l1.p)
    z = l1.direction[3]/l1.direction[1]*(x-getx(l1.p)) + getz(l1.p)
    Point(x, y, z)
    #=!! need exception handlings for cases where direction is 0 !!=#
end

P1=Point(-3.530798,7.196633,0)
dir1=[0.923880,0.382683,0]
P2=Point(-7.196633,3.530798,0)
dir2=[0.382683,0.923880,0]
L1=Line(P1,dir1)
L2=Line(P2,dir2)
findjunction(L1,L2)
