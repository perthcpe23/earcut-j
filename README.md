# earcut-j
Java port of https://github.com/mapbox/earcut

# Usage
Add jar to classpath then
```java
Earcut e = new Earcut();

// coordinates are in order x1, y1, x2, y2, x3, y3,...
int[] indices = e.triangulate(new double[] {10, 0, 0, 50, 60, 60, 70, 10}, new int[] {}, 2);
// result = [1, 0, 3, 3, 2, 1]
// which are indices of two triangles
```
