package RayTracer;

import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.awt.image.BufferedImage;
import javax.imageio.ImageIO;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.sound.sampled.Line;

public class main {

    public static void main(String[] args) {
        // #3
        loadObj("teapot.obj");

        // #9
        Point3D cameraPosition = new Point3D(0, 0, 5);
        Point3D screenCenter = new Point3D(0, 0, 0);
        Point3D screenUp = new Point3D(0, 1, 0);
        double screenWidth = 4.0;
        double screenHeight = 3.0;
        int imageWidth = 800;
        int imageHeight = 600;

        // Initialize screen pixels
        double pixelWidth = screenWidth / imageWidth;
        double pixelHeight = screenHeight / imageHeight;
        List<List<Double>> pixelIntensities = new ArrayList<>();
        for (int y = 0; y < imageHeight; y++) {
            List<Double> rowIntensities = new ArrayList<>();
            for (int x = 0; x < imageWidth; x++) {
                // Calculate ray direction for current pixel
                double screenX = (x - imageWidth / 2.0) * pixelWidth;
                double screenY = (y - imageHeight / 2.0) * pixelHeight;
                Point3D screenPoint = new Point3D(screenX, screenY, 0);
                Point3D rayDirection = subtract(screenPoint, cameraPosition).normalize();

                // Trace ray and calculate illumination intensity
                Line ray = new Line(cameraPosition, rayDirection);
                double intensity = traceRay(ray);

                // Store intensity for current pixel
                rowIntensities.add(intensity);
            }
            pixelIntensities.add(rowIntensities);
        }

        // Print out pixel intensities
        for (int y = 0; y < imageHeight; y++) {
            for (int x = 0; x < imageWidth; x++) {
                System.out.print(pixelIntensities.get(y).get(x) + " ");
            }
            System.out.println();
        }

        // #10
        // 10. Output Screen Pixels as an Image
        // Create BufferedImage to store image
        BufferedImage image = new BufferedImage(imageWidth, imageHeight, BufferedImage.TYPE_INT_RGB);

        // Set pixel intensities to image
        for (int y = 0; y < imageHeight; y++) {
            for (int x = 0; x < imageWidth; x++) {
                // Normalize intensity to RGB range (0-255)
                int intensity = (int) (pixelIntensities.get(y).get(x) * 255);

                // Set pixel color
                int rgb = (intensity << 16) | (intensity << 8) | intensity;
                image.setRGB(x, y, rgb);
            }
        }

        // Output image to file
        File outputFile = new File("output.png");
        try {
            ImageIO.write(image, "png", outputFile);
            System.out.println("Image saved successfully.");
        } catch (IOException e) {
            System.err.println("Error saving image: " + e.getMessage());
        }

        // #3
        for (Face face : faces) {
            for (int i = 0; i < face.vertexIndices.length; i++) {
                Vertex vertex = vertices.get(face.vertexIndices[i]);
                System.out.println("Vertex " + i + ": (" + vertex.x + ", " + vertex.y + ", " + vertex.z + ")");
            }
            System.out.println();
        }

        // #4
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

        Point3D normalize() {
            double magnitude = Math.sqrt(x * x + y * y + z * z);
            if (magnitude == 0) {
                return new Point3D(0, 0, 0);
            }
            return new Point3D(x / magnitude, y / magnitude, z / magnitude);
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
//            Point3D p1 = vertices.get(face.vertexIndices[0]);
//            Point3D p2 = vertices.get(face.vertexIndices[1]);
//            Point3D p3 = vertices.get(face.vertexIndices[2]);
            Point3D p1 = convertTo3D(vertices.get(face.vertexIndices[0]));
            Point3D p2 = convertTo3D(vertices.get(face.vertexIndices[1]));
            Point3D p3 = convertTo3D(vertices.get(face.vertexIndices[2]));

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
    public static Point3D computeTriangleNormal(Point3D p1, Point3D p2, Point3D p3) {
        Point3D v1 = subtract(p2, p1);
        Point3D v2 = subtract(p3, p1);

        // Compute the cross product of v1 and v2
        double normalX = v1.y * v2.z - v1.z * v2.y;
        double normalY = v1.z * v2.x - v1.x * v2.z;
        double normalZ = v1.x * v2.y - v1.y * v2.x;

        // Normalize the result to get the normal vector
        double length = Math.sqrt(normalX * normalX + normalY * normalY + normalZ * normalZ);
        Point3D normal = new Point3D(normalX / length, normalY / length, normalZ / length);

        return normal;
    }

    // 7. Compute the Angle between Two Lines in 3D
    public static double computeAngleBetweenLines(Line line1, Line line2) {
        // Compute the dot product of the direction vectors
        double dotProduct = dotProduct(line1.direction, line2.direction);

        // Compute the magnitudes of the direction vectors
        double magnitude1 = magnitude(line1.direction);
        double magnitude2 = magnitude(line2.direction);

        // Compute the cosine of the angle between the lines
        double cosineAngle = dotProduct / (magnitude1 * magnitude2);

        // Compute the angle in radians
        double angleInRadians = Math.acos(cosineAngle);

        // Convert the angle to degrees
        double angleInDegrees = Math.toDegrees(angleInRadians);

        return angleInDegrees;
    }

    private static double magnitude(Point3D vector) {
        return Math.sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);
    }

    // 8. Compute the Illumination for a Triangle in 3D
    // 9. Raytrace an OBJ in 3D to Identify Illumination Intensities over a Screen
    // of Pixels in 2D

    public static double traceRay(Line ray) {
        // Initialize intensity
        double intensity = 0.0;

        // Check intersections with OBJ geometry
        for (Face face : faces) {
            // Convert face vertices to Point3D objects
//            Point3D p1 = vertices.get(face.vertexIndices[0]);
//            Point3D p2 = vertices.get(face.vertexIndices[1]);
//            Point3D p3 = vertices.get(face.vertexIndices[2]);
            Point3D p1 = convertTo3D(vertices.get(face.vertexIndices[0]));
            Point3D p2 = convertTo3D(vertices.get(face.vertexIndices[1]));
            Point3D p3 = convertTo3D(vertices.get(face.vertexIndices[2]));


            // Create triangle from face vertices
            Triangle triangle = new Triangle(p1, p2, p3);

            // Check if ray intersects with the triangle
            if (isIntersecting(ray, triangle)) {
                // For simplicity, assume constant intensity for intersected triangles
                intensity += 1.0;
            }

            // Compute the normal of the triangle //6
            Point3D normal = computeTriangleNormal(p1, p2, p3);

            // Compute the angle between the ray direction and the normal
            Line normalLine = new Line(p1, normal);
            double angle = computeAngleBetweenLines(ray, normalLine);

            // Compute illumination intensity based on angle
            double illumination = Math.cos(Math.toRadians(angle));
            intensity += illumination;
        }

        return intensity;
    }

    private static Point3D convertTo3D(Vertex vertex) {
        return new Point3D(vertex.x, vertex.y, vertex.z);
    }

    // 10. Output Screen Pixels as an Image

}