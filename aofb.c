#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <unistd.h>
#include <fcntl.h>
#include <linux/fb.h>
#include <linux/fs.h>
#include <sys/mman.h>

#define WIDTH        256
#define HEIGHT       256
#define NSUBSAMPLES  1
#define NAO_SAMPLES  4

#define FBDEV "/dev/fb0"
#define ROWS 4
#define COLS 4
#define FRAMES 700
#define DIV_BYTE 8
#define FRAMES 700

/*
unsigned int wang_hash(unsigned int seed) {
  seed = (seed ^ 61) ^ (seed >> 16);
  seed *= 9;
  seed = seed ^ (seed >> 4);
  seed *= 0x27d4eb2d;
  seed = seed ^ (seed >> 15);
  return seed;
}
*/

// Assume 32bit environment.
/*
inline float fast_sign(float f) {
  float r = 1.0f;
  *((int *)(&r)) |= (*(int *)(&f) & 0x80000000); // mask sign bit in f, set it
                                                 // in r if necessary
  return r;
}
*/

inline float fast_fabs(float x) {
  int f = ((*(int *)(&x)) & 0x7FFFFFFF);
  return *(float *)(&f);
}

// from: http://www.devmaster.net/forums/showthread.php?t=5784
inline float fast_cos(float x) {
  // bring into range [-1,1]
  x = (x + (float)M_PI_2) * (float)M_1_PI;
  // wrap around
  volatile float z =
      (x + 25165824.0f); // must be volatile else it will be optimized out
  x -= (z - 25165824.0f);
  // low precision approximation: return 4.0f * (x - x * fabs(x));
  // higher precision path
  float y = x - x * fast_fabs(x);
  // next iteration for higher precision
  const float Q = 3.1f;
  const float P = 3.6f;
  return y * (Q + P * fast_fabs(y));
}

// from: http://www.devmaster.net/forums/showthread.php?t=5784
inline float fast_sin(float x) {
  // bring into range [-1,1]
  x *= (float)M_1_PI;
  // wrap around
  volatile float z =
      (x + 25165824.0f); // must be volatile else it will be optimized out
  x -= (z - 25165824.0f);
  // low precision approximation: return 4.0f * (x - x * fabs(x));
  // higher precision path
  float y = x - x * fast_fabs(x);
  // next iteration for higher precision
  const float Q = 3.1f;
  const float P = 3.6f;
  return y * (Q + P * fast_fabs(y));
}

static inline float fast_inv_sqrt(float a) {
  union fi32_u {
    float f;
    int32_t i;
  };
  union fi32_u x;
  x.f = a;
  x.i = 0x5f3759df - (x.i >> 1);
  x.f = x.f * (1.5f - x.f * x.f * 0.5f * a);
  x.f = x.f * (1.5f - x.f * x.f * 0.5f * a);
  return x.f;
}

static inline float fast_sqrt(float a) { return fast_inv_sqrt(a) * a; }

/*
static inline float fast_inv(float a) {
  float b = fast_inv_sqrt(fast_fabs(a));
  return fast_sign(a) * b * b;
}
*/

static unsigned long seed_x = 123456789;

void setseed(unsigned long val) { seed_x = val; }

// based on xorshift
float frandom() {
  static unsigned long y = 362436069, z = 521288629, w = 88675123;
  unsigned long x = seed_x;
  unsigned long t;
  t = (x ^ (x << 11));
  seed_x = y;
  y = z;
  z = w;
  w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));

  return (w * (1.0f / 4294967296.0f));
}

typedef struct _vec
{
    float  x;
    float  y;
    float  z;
} vec;


typedef struct _Isect
{
    float  t;
    vec    p;
    vec    n;
    int    hit; 
} Isect;

typedef struct _Sphere
{
    vec    center;
    float  radius;

} Sphere;

typedef struct _Plane
{
    vec    p;
    vec    n;

} Plane;

typedef struct _Ray
{
    vec    org;
    vec    dir;
} Ray;

Sphere spheres[3];
Plane  plane;

static float  vdot(vec v0, vec v1)
{
    return v0.x * v1.x + v0.y * v1.y + v0.z * v1.z;
}

static void vcross(vec *c, vec v0, vec v1)
{
    
    c->x = v0.y * v1.z - v0.z * v1.y;
    c->y = v0.z * v1.x - v0.x * v1.z;
    c->z = v0.x * v1.y - v0.y * v1.x;
}

static void vnormalize(vec *c)
{
    float  length = fast_sqrt(vdot((*c), (*c)));

    if (fast_fabs(length) > 1.0e-6f) {
        //float inv_len = fast_inv(length);
        float inv_len = 1.0f / length;
        c->x *= inv_len;
        c->y *= inv_len;
        c->z *= inv_len;
    }
}

void
ray_sphere_intersect(Isect *isect, const Ray *ray, const Sphere *sphere)
{
    vec rs;

    rs.x = ray->org.x - sphere->center.x;
    rs.y = ray->org.y - sphere->center.y;
    rs.z = ray->org.z - sphere->center.z;

    float B = vdot(rs, ray->dir);
    float C = vdot(rs, rs) - sphere->radius * sphere->radius;
    float D = B * B - C;

    if (D > 0.0) {
        float t = -B - fast_sqrt(D);
        
        if ((t > 0.0) && (t < isect->t)) {
            isect->t = t;
            isect->hit = 1;
            
            isect->p.x = ray->org.x + ray->dir.x * t;
            isect->p.y = ray->org.y + ray->dir.y * t;
            isect->p.z = ray->org.z + ray->dir.z * t;

            isect->n.x = isect->p.x - sphere->center.x;
            isect->n.y = isect->p.y - sphere->center.y;
            isect->n.z = isect->p.z - sphere->center.z;

            vnormalize(&(isect->n));
        }
    }
}

void
ray_plane_intersect(Isect *isect, const Ray *ray, const Plane *plane)
{
    float d = -vdot(plane->p, plane->n);
    float v = vdot(ray->dir, plane->n);

    if (fast_fabs(v) < 1.0e-17) return;

    //float t = -(vdot(ray->org, plane->n) + d) * fast_inv(v);
    float t = -(vdot(ray->org, plane->n) + d) / v;

    if ((t > 0.0) && (t < isect->t)) {
        isect->t = t;
        isect->hit = 1;
        
        isect->p.x = ray->org.x + ray->dir.x * t;
        isect->p.y = ray->org.y + ray->dir.y * t;
        isect->p.z = ray->org.z + ray->dir.z * t;

        isect->n = plane->n;
    }
}

void
orthoBasis(vec *basis, vec n)
{
    basis[2] = n;
    basis[1].x = 0.0; basis[1].y = 0.0; basis[1].z = 0.0;

    if ((n.x < 0.6) && (n.x > -0.6)) {
        basis[1].x = 1.0;
    } else if ((n.y < 0.6) && (n.y > -0.6)) {
        basis[1].y = 1.0;
    } else if ((n.z < 0.6) && (n.z > -0.6)) {
        basis[1].z = 1.0;
    } else {
        basis[1].x = 1.0;
    }

    vcross(&basis[0], basis[1], basis[2]);
    vnormalize(&basis[0]);

    vcross(&basis[1], basis[2], basis[0]);
    vnormalize(&basis[1]);
}


void ambient_occlusion(vec *col, const Isect *isect)
{
    int    i, j;
    int    ntheta = NAO_SAMPLES;
    int    nphi   = NAO_SAMPLES;
    float  eps = 0.0001;

    vec p;

    p.x = isect->p.x + eps * isect->n.x;
    p.y = isect->p.y + eps * isect->n.y;
    p.z = isect->p.z + eps * isect->n.z;

    vec basis[3];
    orthoBasis(basis, isect->n);

    float occlusion = 0.0;

    for (j = 0; j < ntheta; j++) {
        for (i = 0; i < nphi; i++) {
            float  theta = fast_sqrt(frandom());
            float  phi   = 2.0 * M_PI * frandom();

            float  x = fast_cos(phi) * theta;
            float  y = fast_sin(phi) * theta;
            float  z = fast_sqrt(1.0 - theta * theta);

            // local -> global
            float  rx = x * basis[0].x + y * basis[1].x + z * basis[2].x;
            float  ry = x * basis[0].y + y * basis[1].y + z * basis[2].y;
            float  rz = x * basis[0].z + y * basis[1].z + z * basis[2].z;

            Ray ray;

            ray.org = p;
            ray.dir.x = rx;
            ray.dir.y = ry;
            ray.dir.z = rz;

            Isect occIsect;
            occIsect.t   = 1.0e+17;
            occIsect.hit = 0;

            ray_sphere_intersect(&occIsect, &ray, &spheres[0]); 
            ray_sphere_intersect(&occIsect, &ray, &spheres[1]); 
            ray_sphere_intersect(&occIsect, &ray, &spheres[2]); 
            ray_plane_intersect (&occIsect, &ray, &plane); 

            if (occIsect.hit) occlusion += 1.0;
            
        }
    }

    occlusion = (ntheta * nphi - occlusion) / (float)(ntheta * nphi);

    col->x = occlusion;
    col->y = occlusion;
    col->z = occlusion;
}

unsigned char
clamp(float  f)
{
  int i = (int)(f * 255.5);

  if (i < 0) i = 0;
  if (i > 255) i = 255;

  return (unsigned char)i;
}


#if 0
void
render(unsigned char *img, int w, int h, int nsubsamples)
{
    int x, y;
    int u, v;

    float *fimg = (float  *)malloc(sizeof(float ) * w * h * 3);
    memset((void *)fimg, 0, sizeof(float ) * w * h * 3);

    for (y = 0; y < h; y++) {
        for (x = 0; x < w; x++) {
            
            for (v = 0; v < nsubsamples; v++) {
                for (u = 0; u < nsubsamples; u++) {
                    float px = (x + (u / (float)nsubsamples) - (w / 2.0)) / (w / 2.0);
                    float py = -(y + (v / (float)nsubsamples) - (h / 2.0)) / (h / 2.0);

                    Ray ray;

                    ray.org.x = 0.0;
                    ray.org.y = 0.0;
                    ray.org.z = 0.0;

                    ray.dir.x = px;
                    ray.dir.y = py;
                    ray.dir.z = -1.0;
                    vnormalize(&(ray.dir));

                    Isect isect;
                    isect.t   = 1.0e+17;
                    isect.hit = 0;

                    ray_sphere_intersect(&isect, &ray, &spheres[0]);
                    ray_sphere_intersect(&isect, &ray, &spheres[1]);
                    ray_sphere_intersect(&isect, &ray, &spheres[2]);
                    ray_plane_intersect (&isect, &ray, &plane);

                    if (isect.hit) {
                        vec col;
                        ambient_occlusion(&col, &isect);

                        fimg[3 * (y * w + x) + 0] += col.x;
                        fimg[3 * (y * w + x) + 1] += col.y;
                        fimg[3 * (y * w + x) + 2] += col.z;
                    }

                }
            }

            fimg[3 * (y * w + x) + 0] /= (float)(nsubsamples * nsubsamples);
            fimg[3 * (y * w + x) + 1] /= (float)(nsubsamples * nsubsamples);
            fimg[3 * (y * w + x) + 2] /= (float)(nsubsamples * nsubsamples);
        
            img[3 * (y * w + x) + 0] = clamp(fimg[3 *(y * w + x) + 0]);
            img[3 * (y * w + x) + 1] = clamp(fimg[3 *(y * w + x) + 1]);
            img[3 * (y * w + x) + 2] = clamp(fimg[3 *(y * w + x) + 2]);
        }
    }

}
#endif

void
init_scene(float timeval)
{
#if 0
    spheres[0].center.x = -2.0;
    spheres[0].center.y =  0.0;
    spheres[0].center.z = -3.5;
    spheres[0].radius = 0.5;
    
    spheres[1].center.x = -0.5;
    spheres[1].center.y =  0.0;
    spheres[1].center.z = -3.0;
    spheres[1].radius = 0.5;
    
    spheres[2].center.x =  1.0;
    spheres[2].center.y =  0.0;
    spheres[2].center.z = -2.2;
    spheres[2].radius = 0.5;
#else
  spheres[0].center.x = -1.0f + fast_sin(timeval);
  spheres[0].center.y = 1.0f;
  spheres[0].center.z = -4.5f + fast_cos(timeval);
  spheres[0].radius = 1.5f;

  spheres[1].center.x = 1.5f + fast_cos(timeval);
  spheres[1].center.y = 0.5f;
  spheres[1].center.z = -3.5f + fast_sin(timeval);
  spheres[1].radius = 1.0f;

  spheres[2].center.x = 0.5f + fast_cos(timeval * 2.3f);
  spheres[2].center.y = 0.0f;
  spheres[2].center.z = -3.2f + fast_sin(timeval * 2.3f);
  spheres[2].radius = 0.5f;

#endif


    plane.p.x = 0.0;
    plane.p.y = -0.5;
    plane.p.z = 0.0;

    plane.n.x = 0.0;
    plane.n.y = 1.0;
    plane.n.z = 0.0;

}

/*
void
saveppm(const char *fname, int w, int h, unsigned char *img)
{
    FILE *fp;

    fp = fopen(fname, "wb");
    assert(fp);

    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", w, h);
    fprintf(fp, "255\n");
    fwrite(img, w * h * 3, 1, fp);
    fclose(fp);
}
*/

int
main(int argc, char **argv)
{
    unsigned int *img;

    int fb;
    char *fbptr;
    struct fb_var_screeninfo vinfo;
    struct fb_fix_screeninfo finfo;
    unsigned int frame = 1;

    unsigned int w = WIDTH;
    unsigned int h = HEIGHT;
    int xres,yres,vbpp,line_len;
    size_t screensize;
    long int location;
    int ret;

    if( argc <= 1 ) {
        fprintf(stderr, "usage: %s /dev/fb[0-9]\n", argv[0]);
        exit(1);
    }

    img = (unsigned int *)malloc(WIDTH * HEIGHT * 4);
    if( !img ) {
        fprintf(stderr, "malloc error\n");
        exit(1);
    }

    /* fb = open(FBDEV, O_RDWR); */
    fb = open(argv[1], O_RDWR);
    if( !fb ) {
        fprintf(stderr, "fb open error\n");
        exit(1);
    }

    if((ret = ioctl(fb, FBIOGET_FSCREENINFO , &finfo )) != 0) {
        fprintf(stderr, "fixed information error!: %d\n", ret);
        exit(2);
    }
    printf("smem_start: %lx\n", finfo.smem_start);
    printf("smem_len: %d\n", finfo.smem_len);
    printf("line_length: %d\n", finfo.line_length);

    if ((ret = ioctl( fb, FBIOGET_VSCREENINFO , &vinfo )) != 0) {
        fprintf(stderr, "variable information error!: %d\n", ret);
        exit(3);
    }
    printf("xres: %d\n", vinfo.xres);
    printf("yres: %d\n", vinfo.yres);
    printf("xres_virtual: %d\n", vinfo.xres_virtual);
    printf("yres_virtual: %d\n", vinfo.yres_virtual);
    printf("xoffset: %d\n", vinfo.xoffset);
    printf("yoffset: %d\n", vinfo.yoffset);
    printf("bit per pixel: %d\n", vinfo.bits_per_pixel);

    xres = vinfo.xres;
    yres = vinfo.yres;
    vbpp = vinfo.bits_per_pixel;
    line_len = finfo.line_length;

    /* screen size in byte */
    screensize = xres * yres * vbpp / DIV_BYTE;

    fbptr = (char *)mmap(0, screensize, PROT_READ | PROT_WRITE,MAP_SHARED, fb, 0);
    if( (int)fbptr == -1 ){
        fprintf(stderr, "cannot mmap framebuffer\n");
        exit(4);
    }

    while(1){
        int x, y;
        frame++;
        init_scene((float)frame * 0.1f);

        /* display */
        for(y = 0; y < HEIGHT; y++) {
            for(x = 0; x < WIDTH; x++) {
                int u, v;
                float fimg[3] ={ 0.0f, 0.0f, 0.0f };

                for (v = 0; v < NSUBSAMPLES; v++) {
                    for (u = 0; u < NSUBSAMPLES; u++) {
                        float px = (x + (u / (float)NSUBSAMPLES) - (w / 2.0)) / (w / 2.0);
                        float py = -(y + (v / (float)NSUBSAMPLES) - (h / 2.0)) / (h / 2.0);

                        Ray ray;

                        ray.org.x = 0.0f;
                        ray.org.y = 0.0f;
                        ray.org.z = 1.0f;

                        ray.dir.x = px;
                        ray.dir.y = py;
                        ray.dir.z = -1.0f;
                        vnormalize(&(ray.dir));

                        Isect isect;
                        isect.t   = 1.0e+6f;
                        isect.hit = 0;

                        ray_sphere_intersect(&isect, &ray, &spheres[0]);
                        ray_sphere_intersect(&isect, &ray, &spheres[1]);
                        ray_sphere_intersect(&isect, &ray, &spheres[2]);
                        ray_plane_intersect (&isect, &ray, &plane);

                        if (isect.hit) {
                            vec col;
                            ambient_occlusion(&col, &isect);

                            fimg[0] += col.x;
                            fimg[1] += col.y;
                            fimg[2] += col.z;
                        }

                    }
                }

                location = ((128 + x + vinfo.xoffset) * vbpp / DIV_BYTE) +
                    (128 + y + vinfo.yoffset) * line_len;

                fimg[0] /= (float)(NSUBSAMPLES * NSUBSAMPLES);
                fimg[1] /= (float)(NSUBSAMPLES * NSUBSAMPLES);
                fimg[2] /= (float)(NSUBSAMPLES * NSUBSAMPLES);

                /* ARGB */
                *((unsigned int*)(fbptr + location)) = (clamp(fimg[0]) << 16) |
                                                       (clamp(fimg[1]) << 8) |
                                                       (clamp(fimg[2]));
            }
        }
        printf(".");
        fflush(stdout);
        if (frame > FRAMES) {
            break;
        }
    }

    munmap(fbptr , screensize);
    close(fb);

    return 0;
}
