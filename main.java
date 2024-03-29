package RayTracer;

import java.io.File;
import java.io.IOException;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.sound.sampled.Line;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class Main {

    public static void main(String[] args) {
        loadObj("teapot.obj");

        for (Face face : faces) {
            for (int i = 0; i < face.vertexIndices.length; i++) {
                Vertex vertex = vertices.get(face.vertexIndices[i]);
                System.out.println("Vertex " + i + ": (" + vertex.x + ", " + vertex.y + ", " + vertex.z + ")");
            }
            System.out.println();
        }

        Point3D p1 = new Point3D(0, 0, 0);
        Point3D p2 = new Point3D(1, 0, 0);
        Point3D p3 = new Point3D(0, 1, 0);
        Triangle triangle = new Triangle(p1, p2, p3);

        Point3D origin = new Point3D(0.5, 0.5, 1);
        Point3D direction = new Point3D(0, 0, -1);
        Line line = new Line(origin, direction);

        boolean isIntersecting = isIntersecting(line, triangle);
        System.out.println("Is intersecting: " + isIntersecting);

        List<Face> intersectingFaces = intersectingFaces(line);
        System.out.println("Intersecting Faces: " + intersectingFaces);
    }

    static List<Vertex> vertices = new ArrayList<>();
    static List<Face> faces = new ArrayList<>();

    public static void loadObj(String filename) {
        try (BufferedReader br = new BufferedReader(new FileReader(filename))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("v ")) {
                    String[] parts = line.split("\\s+");
                    float x = Float.parseFloat(parts[1]);
                    float y = Float.parseFloat(parts[2]);
                    float z = Float.parseFloat(parts[3]);
                    vertices.add(new Vertex(x, y, z));
                } else if (line.startsWith("f ")) {
                    String[] parts = line.split("\\s+");
                    int[] vertexIndices = new int[parts.length - 1];
                    for (int i = 1; i < parts.length; i++) {
                        String[] indices = parts[i].split("/");
                        int vertexIndex = Integer.parseInt(indices[0]) - 1;
                        vertexIndices[i - 1] = vertexIndex;
                    }
                    faces.add(new Face(vertexIndices));
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    static class Vertex {
        public float x, y, z;

        public Vertex(float x, float y, float z) {
            this.x = x;
            this.y = y;
            this.z = z;
        }
    }

    static class Face {
        public int[] vertexIndices;

        public Face(int[] vertexIndices) {
            this.vertexIndices = vertexIndices;
        }
    }

    static class Point3D {
        double x, y, z;

        public Point3D(double x, double y, double z) {
            this.x = x;
            this.y = y;
            this.z = z;
        }
    }

    static class Triangle {
        Point3D p1, p2, p3;

        public Triangle(Point3D p1, Point3D p2, Point3D p3) {
            this.p1 = p1;
            this.p2 = p2;
            this.p3 = p3;
        }
    }

    static class Line {
        Point3D origin, direction;

        public Line(Point3D origin, Point3D direction) {
            this.origin = origin;
            this.direction = direction;
        }
    }

    // 4. Determine whether a Line and a Triangle Intersect in 3D

    public static boolean isIntersecting(Line line, Triangle triangle) {
        Point3D e1 = subtract(triangle.p2, triangle.p1);
        Point3D e2 = subtract(triangle.p3, triangle.p1);
        Point3D h = crossProduct(line.direction, e2);
        double a = dotProduct(e1, h);

        if (a > -0.00001 && a < 0.00001)
            return false;

        double f = 1.0 / a;
        Point3D s = subtract(line.origin, triangle.p1);
        double u = f * dotProduct(s, h);

        if (u < 0.0 || u > 1.0)
            return false;

        Point3D q = crossProduct(s, e1);
        double v = f * dotProduct(line.direction, q);

        if (v < 0.0 || u + v > 1.0)
            return false;

        double t = f * dotProduct(e2, q);

        return t > 0.00001;
    }

    private static Point3D subtract(Point3D a, Point3D b) {
        return new Point3D(a.x - b.x, a.y - b.y, a.z - b.z);
    }

    private static double dotProduct(Point3D a, Point3D b) {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    private static Point3D crossProduct(Point3D a, Point3D b) {
        double x = a.y * b.z - a.z * b.y;
        double y = a.z * b.x - a.x * b.z;
        double z = a.x * b.y - a.y * b.x;
        return new Point3D(x, y, z);
    }

    // 5. Determine Which Faces from an OBJ File a Given Line Intersects With in 3D

    public static List<Face> intersectingFaces(Line line) {
        List<Face> intersecting = new ArrayList<>();
        for (Face face : faces) {
            // Assume each face is a triangle
            Point3D p1 = vertices.get(face.vertexIndices[0]);
            Point3D p2 = vertices.get(face.vertexIndices[1]);
            Point3D p3 = vertices.get(face.vertexIndices[2]);
            Triangle triangle = new Triangle(p1, p2, p3);
            if (isIntersecting(line, triangle)) {
                intersecting.add(face);
            }
        }
        return intersecting;
    }

    // 2. Compute Absolute Value

    public static Double ABSValue(Double input) {
        Double ret = Math.abs(input);
        return ret;
    }

    // 3. Compute the Distance between a Point and a Line in 3D

    public static float distancePointLine(Vertex point, Vertex linePoint1, Vertex linePoint2) {
        float x0 = point.x;
        float y0 = point.y;
        float z0 = point.z;

        float x1 = linePoint1.x;
        float y1 = linePoint1.y;
        float z1 = linePoint1.z;

        float x2 = linePoint2.x;
        float y2 = linePoint2.y;
        float z2 = linePoint2.z;

        float dx = x2 - x1;
        float dy = y2 - y1;
        float dz = z2 - z1;

        float t = ((x0 - x1) * dx + (y0 - y1) * dy + (z0 - z1) * dz) / (dx * dx + dy * dy + dz * dz);

        float closestX, closestY, closestZ;

        if (t < 0) {
            closestX = x1;
            closestY = y1;
            closestZ = z1;
        } else if (t > 1) {
            closestX = x2;
            closestY = y2;
            closestZ = z2;
        } else {
            closestX = x1 + t * dx;
            closestY = y1 + t * dy;
            closestZ = z1 + t * dz;
        }

        float distance = (float) Math.sqrt((x0 - closestX) * (x0 - closestX) +
                (y0 - closestY) * (y0 - closestY) +
                (z0 - closestZ) * (z0 - closestZ));
        return distance;
    }

    // 6. Compute the Normal of a Triangle in 3D
    // 7. Compute the Angle between Two Lines in 3D




    // 8. Compute the Illumination for a Triangle in 3D
    static class Triangle {
    Point3D p1, p2, p3;
    Vector3D normal; // Normal vector of the triangle

    public Triangle(Point3D p1, Point3D p2, Point3D p3) {
        this.p1 = p1;
        this.p2 = p2;
        this.p3 = p3;
        this.normal = calculateNormal();
    }

    // Calculate the normal vector of the triangle
    private Vector3D calculateNormal() {
        Vector3D e1 = new Vector3D(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);
        Vector3D e2 = new Vector3D(p3.x - p1.x, p3.y - p1.y, p3.z - p1.z);
        return e1.crossProduct(e2).normalize();
    }

    // Compute the illumination intensity at a given point on the triangle
    public double computeIllumination(Vector3D lightDirection, Vector3D viewDirection, double ambientIntensity, 
                                      double diffuseIntensity, double specularIntensity, double shininess) {
        // Ambient component
        double ambient = ambientIntensity;

        // Diffuse component
        double diffuse = Math.max(0, normal.dotProduct(lightDirection));

        // Specular component
        Vector3D reflection = normal.scale(2 * normal.dotProduct(lightDirection)).subtract(lightDirection).normalize();
        double specular = Math.pow(Math.max(0, reflection.dotProduct(viewDirection)), shininess);

        // Total intensity
        return ambient + diffuseIntensity * diffuse + specularIntensity * specular;
    }
}

    // 9. Raytrace an OBJ in 3D to Identify Illumination Intensities over a Screen of Pixels in 2D



    // 10. Output Screen Pixels as an Image

}