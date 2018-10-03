package com.perth.mapbox.earcut;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

class Node{
	// vertex index in coordinates array
	int i;
	
	// vertex coordinates
	double x, y;
	
	// z-order curve value
	double z = Double.MIN_VALUE;
	
	// indicates whether this is a steiner point
	boolean steiner = false;
	
	// previous and next vertex nodes in a polygon ring
	Node prev, next;
	
	// previous and next nodes in z-order
	Node prevZ, nextZ;
	
	Node(int i, double x, double y){
		this.i = i;
		this.x = x;
		this.y = y;
	}
}

public class Earcut {

	private boolean hasHoles;
	private int outerLen;
	private Node outerNode;
	private List<Integer> triangles;
	        
	public Earcut() {}
	
	public int[] triangulate(double[] data, int[] holeIndices, int dim){
		hasHoles = holeIndices != null && holeIndices.length > 0;
		outerLen = hasHoles ? holeIndices.length * dim : data.length;
		outerNode = linkedList(data, 0, outerLen, dim, true);
		triangles = new ArrayList<>();
		
		if (outerNode == null) return new int[] {};

	    double minX = 0, minY = 0, maxX, maxY, x, y;
	    int invSize = 0;

	    if (hasHoles) outerNode = eliminateHoles(data, holeIndices, outerNode, dim);

	    // if the shape is not too simple, we'll use z-order curve hash later; calculate polygon bbox
	    if (data.length > 80 * dim) {
	        minX = maxX = data[0];
	        minY = maxY = data[1];

	        for (int i = dim; i < outerLen; i += dim) {
	            x = data[i];
	            y = data[i + 1];
	            if (x < minX) minX = x;
	            if (y < minY) minY = y;
	            if (x > maxX) maxX = x;
	            if (y > maxY) maxY = y;
	        }

 	        // minX, minY and invSize are later used to transform coords into integers for z-order calculation
	        invSize = (int)Math.max(maxX - minX, maxY - minY);
	        invSize = invSize != 0 ? 1 / invSize : 0;
	    }

	    earcutLinked(outerNode, triangles, dim, minX, minY, invSize, 0);
	    
	    int[] indices = new int[triangles.size()];
	    for(int i=0;i<indices.length;i++) {
	    	indices[i] = triangles.get(i);
	    }
	    
	    return indices;
	}
	
	Node linkedList(double[] data, int start, int end, int dim, boolean clockwise) {
		int i;
		Node last = null;

	    if (clockwise == (signedArea(data, start, end, dim) > 0)) {
	        for (i = start; i < end; i += dim) last = insertNode(i, data[i], data[i + 1], last);
	    } else {
	        for (i = end - dim; i >= start; i -= dim) last = insertNode(i, data[i], data[i + 1], last);
	    }

	    if (last != null && equals(last, last.next)) {
	        removeNode(last);
	        last = last.next;
	    }

	    return last;
	}
	
	void removeNode(Node p) {
		p.next.prev = p.prev;
	    p.prev.next = p.next;

	    if (p.prevZ != null) p.prevZ.nextZ = p.nextZ;
	    if (p.nextZ != null) p.nextZ.prevZ = p.prevZ;
	}
	
	// create a node and optionally link it with previous one (in a circular doubly linked list)
	Node insertNode(int i, double x, double y, Node last) {
	    Node p = new Node(i, x, y);

	    if (last == null) {
	        p.prev = p;
	        p.next = p;

	    } else {
	        p.next = last.next;
	        p.prev = last;
	        last.next.prev = p;
	        last.next = p;
	    }
	    
	    return p;
	}

	// eliminate colinear or duplicate points
	Node filterPoints(Node start, Node end) {
	    if (start == null) return start;
	    if (end == null) end = start;

	    Node p = start;
	    boolean again;
	    do {
	        again = false;

	        if (!p.steiner && (equals(p, p.next) || area(p.prev, p, p.next) == 0)) {
	            removeNode(p);
	            p = end = p.prev;
	            if (p == p.next) break;
	            again = true;

	        } else {
	            p = p.next;
	        }
	    } while (again || p != end);

	    return end;
	}

	// main ear slicing loop which triangulates a polygon (given as a linked list)
	void earcutLinked(Node ear, List<Integer> triangles, int dim, double minX, double minY, int invSize, int pass) {
	    if (ear == null) return;

	    // interlink polygon nodes in z-order
	    if (pass == 0 && invSize != 0) indexCurve(ear, minX, minY, invSize);

	    Node stop = ear, prev, next;

	    // iterate through ears, slicing them one by one
	    while (ear.prev != ear.next) {
	        prev = ear.prev;
	        next = ear.next;

	        if (invSize != 0 ? isEarHashed(ear, minX, minY, invSize) : isEar(ear)) {
	            // cut off the triangle
	            triangles.add(prev.i / dim);
	            triangles.add(ear.i / dim);
	            triangles.add(next.i / dim);

	            removeNode(ear);

	            // skipping the next vertex leads to less sliver triangles
	            ear = next.next;
	            stop = next.next;

	            continue;
	        }

	        ear = next;

	        // if we looped through the whole remaining polygon and can't find any more ears
	        if (ear == stop) {
	            // try filtering points and slicing again
	            if (pass == 0) {
	                earcutLinked(filterPoints(ear,null), triangles, dim, minX, minY, invSize, 1);

	            // if this didn't work, try curing all small self-intersections locally
	            } else if (pass == 1) {
	                ear = cureLocalIntersections(ear, triangles, dim);
	                earcutLinked(ear, triangles, dim, minX, minY, invSize, 2);

	            // as a last resort, try splitting the remaining polygon into two
	            } else if (pass == 2) {
	                splitEarcut(ear, triangles, dim, minX, minY, invSize);
	            }

	            break;
	        }
	    }
	}
	
	// check whether a polygon node forms a valid ear with adjacent nodes
	boolean isEar(Node ear) {
	    Node a = ear.prev, b = ear, c = ear.next;

	    if (area(a, b, c) >= 0) return false; // reflex, can't be an ear

	    // now make sure we don't have other points inside the potential ear
	    Node p = ear.next.next;

	    while (p != ear.prev) {
	        if (pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y) &&
	            area(p.prev, p, p.next) >= 0) return false;
	        p = p.next;
	    }

	    return true;
	}

	boolean isEarHashed(Node ear, double minX, double minY, int invSize) {
	    Node a = ear.prev, b = ear, c = ear.next;

	    if (area(a, b, c) >= 0) return false; // reflex, can't be an ear

	    // triangle bbox; min & max are calculated like this for speed
	    double minTX = a.x < b.x ? (a.x < c.x ? a.x : c.x) : (b.x < c.x ? b.x : c.x);
	    double minTY = a.y < b.y ? (a.y < c.y ? a.y : c.y) : (b.y < c.y ? b.y : c.y);
	    double maxTX = a.x > b.x ? (a.x > c.x ? a.x : c.x) : (b.x > c.x ? b.x : c.x);
	    double maxTY = a.y > b.y ? (a.y > c.y ? a.y : c.y) : (b.y > c.y ? b.y : c.y);

	    // z-order range for the current triangle bbox;
	    double minZ = zOrder((int)minTX, (int)minTY, minX, minY, invSize);
	    double maxZ = zOrder((int)maxTX, (int)maxTY, minX, minY, invSize);

	    Node p = ear.prevZ, n = ear.nextZ;

	    // look for points inside the triangle in both directions
	    while (p != null && p.z >= minZ && n != null && n.z <= maxZ) {
	        if (p != ear.prev && p != ear.next &&
	            pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y) &&
	            area(p.prev, p, p.next) >= 0) return false;
	        p = p.prevZ;

	        if (n != ear.prev && n != ear.next &&
	            pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, n.x, n.y) &&
	            area(n.prev, n, n.next) >= 0) return false;
	        n = n.nextZ;
	    }

	    // look for remaining points in decreasing z-order
	    while (p != null && p.z >= minZ) {
	        if (p != ear.prev && p != ear.next &&
	            pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y) &&
	            area(p.prev, p, p.next) >= 0) return false;
	        p = p.prevZ;
	    }

	    // look for remaining points in increasing z-order
	    while (n != null && n.z <= maxZ) {
	        if (n != ear.prev && n != ear.next &&
	            pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, n.x, n.y) &&
	            area(n.prev, n, n.next) >= 0) return false;
	        n = n.nextZ;
	    }

	    return true;
	}

	// go through all polygon nodes and cure small local self-intersections
	Node cureLocalIntersections(Node start, List<Integer> triangles, int dim) {
	    Node p = start;
	    do {
	        Node a = p.prev;
    		Node b = p.next.next;

	        if (!equals(a, b) && intersects(a, p, p.next, b) && locallyInside(a, b) && locallyInside(b, a)) {

	            triangles.add(a.i / dim);
	            triangles.add(p.i / dim);
	            triangles.add(b.i / dim);

	            // remove two nodes involved
	            removeNode(p);
	            removeNode(p.next);

	            p = start = b;
	        }
	        p = p.next;
	    } while (p != start);

	    return p;
	}

	// try splitting polygon into two and triangulate them independently
	void splitEarcut(Node start, List<Integer> triangles, int dim, double minX, double minY, int invSize) {
	    // look for a valid diagonal that divides the polygon into two
	    Node a = start;
	    do {
	        Node b = a.next.next;
	        while (b != a.prev) {
	            if (a.i != b.i && isValidDiagonal(a, b)) {
	                // split the polygon in two by the diagonal
	                Node c = splitPolygon(a, b);

	                // filter colinear points around the cuts
	                a = filterPoints(a, a.next);
	                c = filterPoints(c, c.next);

	                // run earcut on each half
	                earcutLinked(a, triangles, dim, minX, minY, invSize, 0);
	                earcutLinked(c, triangles, dim, minX, minY, invSize, 0);
	                return;
	            }
	            b = b.next;
	        }
	        a = a.next;
	    } while (a != start);
	}

	// link every hole into the outer loop, producing a single-ring polygon without holes
	Node eliminateHoles(double[] data, int[] holeIndices, Node outerNode, int dim) {
	    List<Node> queue = new ArrayList<>();
	    int i, len, start, end;
	    Node list;

	    for (i = 0, len = holeIndices.length; i < len; i++) {
	        start = holeIndices[i] * dim;
	        end = i < len - 1 ? holeIndices[i + 1] * dim : data.length;
	        list = linkedList(data, start, end, dim, false);
	        if (list == list.next) list.steiner = true;
	        queue.add(getLeftmost(list));
	    }

	    Collections.sort(queue,new Comparator<Node>() {

			@Override
			public int compare(Node a, Node b) {
				return compareX(a, b);
			}
	    	
		});

	    // process holes from left to right
	    for (i = 0; i < queue.size(); i++) {
	        eliminateHole(queue.get(i), outerNode);
	        outerNode = filterPoints(outerNode, outerNode.next);
	    }

	    return outerNode;
	}

	int compareX(Node a, Node b) {
	    return (int) (a.x - b.x);
	}

	// find a bridge between vertices that connects hole with an outer ring and and link it
	void eliminateHole(Node hole, Node outerNode) {
	    outerNode = findHoleBridge(hole, outerNode);
	    if (outerNode != null) {
	        Node b = splitPolygon(outerNode, hole);
	        filterPoints(b, b.next);
	    }
	}
	
	// David Eberly's algorithm for finding a bridge between hole and outer polygon
	Node findHoleBridge(Node hole, Node outerNode) {
	    Node p = outerNode;
	    double hx = hole.x;
	    double hy = hole.y;
	    double qx = Double.MIN_VALUE;
	    Node m = null;

	    // find a segment intersected by a ray from the hole's leftmost point to the left;
	    // segment's endpoint with lesser x will be potential connection point
	    do {
	        if (hy <= p.y && hy >= p.next.y && p.next.y != p.y) {
	            double x = p.x + (hy - p.y) * (p.next.x - p.x) / (p.next.y - p.y);
	            if (x <= hx && x > qx) {
	                qx = x;
	                if (x == hx) {
	                    if (hy == p.y) return p;
	                    if (hy == p.next.y) return p.next;
	                }
	                m = p.x < p.next.x ? p : p.next;
	            }
	        }
	        p = p.next;
	    } while (p != outerNode);

	    if (m == null) return null;

	    if (hx == qx) return m.prev; // hole touches outer segment; pick lower endpoint

	    // look for points inside the triangle of hole point, segment intersection and endpoint;
	    // if there are no points found, we have a valid connection;
	    // otherwise choose the point of the minimum angle with the ray as connection point

	    Node stop = m;
	    double mx = m.x,
	        my = m.y,
	        tanMin = Double.MAX_VALUE,
	        tan;

	    p = m.next;

	    while (p != stop) {
	        if (hx >= p.x && p.x >= mx && hx != p.x &&
	                pointInTriangle(hy < my ? hx : qx, hy, mx, my, hy < my ? qx : hx, hy, p.x, p.y)) {

	            tan = Math.abs(hy - p.y) / (hx - p.x); // tangential

	            if ((tan < tanMin || (tan == tanMin && p.x > m.x)) && locallyInside(p, hole)) {
	                m = p;
	                tanMin = tan;
	            }
	        }

	        p = p.next;
	    }

	    return m;
	}

	// link two polygon vertices with a bridge; if the vertices belong to the same ring, it splits polygon into two;
	// if one belongs to the outer ring and another to a hole, it merges it into a single ring
	Node splitPolygon(Node a, Node b) {
	    Node a2 = new Node(a.i, a.x, a.y);
	    Node b2 = new Node(b.i, b.x, b.y);
	    Node an = a.next;
	    Node bp = b.prev;

	    a.next = b;
	    b.prev = a;

	    a2.next = an;
	    an.prev = a2;

	    b2.next = a2;
	    a2.prev = b2;

	    bp.next = b2;
	    b2.prev = bp;

	    return b2;
	}
	
	// check if the middle point of a polygon diagonal is inside the polygon
	boolean middleInside(Node a, Node b) {
	    Node p = a;
	    boolean inside = false;
	    double px = (a.x + b.x) / 2;
	    double py = (a.y + b.y) / 2;
	    
	    do {
	        if (((p.y > py) != (p.next.y > py)) && p.next.y != p.y &&
	                (px < (p.next.x - p.x) * (py - p.y) / (p.next.y - p.y) + p.x))
	            inside = !inside;
	        p = p.next;
	    } while (p != a);

	    return inside;
	}
		
	// check if a polygon diagonal is locally inside the polygon
	boolean locallyInside(Node a, Node b) {
	    return area(a.prev, a, a.next) < 0 ?
	        area(a, b, a.next) >= 0 && area(a, a.prev, b) >= 0 :
	        area(a, b, a.prev) < 0 || area(a, a.next, b) < 0;
	}
	
	// check if a polygon diagonal intersects any polygon segments
	boolean intersectsPolygon(Node a, Node b) {
	    Node p = a;
	    do {
	        if (p.i != a.i && p.next.i != a.i && p.i != b.i && p.next.i != b.i && intersects(p, p.next, a, b)) 
	        	return true;
	        p = p.next;
	    } while (p != a);

	    return false;
	}
	
	// check if two segments intersect
	boolean intersects(Node p1, Node q1, Node p2, Node q2) {
	    if ((equals(p1, q1) && equals(p2, q2)) ||
	        (equals(p1, q2) && equals(p2, q1))) return true;
	    
	    return area(p1, q1, p2) > 0 != area(p1, q1, q2) > 0 &&
	           area(p2, q2, p1) > 0 != area(p2, q2, q1) > 0;
	}
	
	// z-order of a point given coords and inverse of the longer side of data bbox
	int zOrder(int x, int y, double minX, double minY, int invSize) {
	    // coords are transformed into non-negative 15-bit integer range
	    x = (int) (32767 * (x - minX) * invSize);
	    y = (int) (32767 * (y - minY) * invSize);

	    x = (x | (x << 8)) & 0x00FF00FF;
	    x = (x | (x << 4)) & 0x0F0F0F0F;
	    x = (x | (x << 2)) & 0x33333333;
	    x = (x | (x << 1)) & 0x55555555;

	    y = (y | (y << 8)) & 0x00FF00FF;
	    y = (y | (y << 4)) & 0x0F0F0F0F;
	    y = (y | (y << 2)) & 0x33333333;
	    y = (y | (y << 1)) & 0x55555555;

	    return x | (y << 1);
	}

	// find the leftmost node of a polygon ring
	Node getLeftmost(Node start) {
	    Node p = start, leftmost = start;
	    do {
	        if (p.x < leftmost.x) leftmost = p;
	        p = p.next;
	    } while (p != start);

	    return leftmost;
	}

	// check if a point lies within a convex triangle
	boolean pointInTriangle(double ax, double ay, double bx, double by, double cx, double cy, double px, double py) {
	    return (cx - px) * (ay - py) - (ax - px) * (cy - py) >= 0 &&
	           (ax - px) * (by - py) - (bx - px) * (ay - py) >= 0 &&
	           (bx - px) * (cy - py) - (cx - px) * (by - py) >= 0;
	}

	// check if a diagonal between two polygon nodes is valid (lies in polygon interior)
	boolean isValidDiagonal(Node a, Node b) {
	    return a.next.i != b.i && a.prev.i != b.i && !intersectsPolygon(a, b) &&
	           locallyInside(a, b) && locallyInside(b, a) && middleInside(a, b);
	}

	// signed area of a triangle
	double area(Node p, Node q, Node r) {
	    return (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
	}

	// check if two points are equal
	boolean equals(Node p1, Node p2) {
	    return p1.x == p2.x && p1.y == p2.y;
	}
	
	// interlink polygon nodes in z-order
	void indexCurve(Node start, double minX, double minY, int invSize) {
	    Node p = start;
	    do {
	        if (p.z == Double.MIN_VALUE) p.z = zOrder((int)p.x, (int)p.y, minX, minY, invSize);
	        p.prevZ = p.prev;
	        p.nextZ = p.next;
	        p = p.next;
	    } while (p != start);

	    p.prevZ.nextZ = null;
	    p.prevZ = null;

	    sortLinked(p);
	}

	// Simon Tatham's linked list merge sort algorithm
	// http://www.chiark.greenend.org.uk/~sgtatham/algorithms/listsort.html
	Node sortLinked(Node list) {
	    int i;
	    Node p, q, e, tail;
	    int numMerges, pSize, qSize;
	    int inSize = 1;

	    do {
	        p = list;
	        list = null;
	        tail = null;
	        numMerges = 0;

	        while (p != null) {
	            numMerges++;
	            q = p;
	            pSize = 0;
	            for (i = 0; i < inSize; i++) {
	                pSize++;
	                q = q.nextZ;
	                if (q == null) break;
	            }
	            qSize = inSize;

	            while (pSize > 0 || (qSize > 0 && q != null)) {

	                if (pSize != 0 && (qSize == 0 || q == null || p.z <= q.z)) {
	                    e = p;
	                    p = p.nextZ;
	                    pSize--;
	                } else {
	                    e = q;
	                    q = q.nextZ;
	                    qSize--;
	                }

	                if (tail != null) tail.nextZ = e;
	                else list = e;

	                e.prevZ = tail;
	                tail = e;
	            }

	            p = q;
	        }

	        tail.nextZ = null;
	        inSize *= 2;

	    } while (numMerges > 1);

	    return list;
	}
	
	double signedArea(double[] data, int start, int end, int dim) {
	    double sum = 0;
	    for (int i = start, j = end - dim; i < end; i += dim) {
	        sum += (data[j] - data[i]) * (data[i + 1] + data[j + 1]);
	        j = i;
	    }
	    return sum;
	}
}