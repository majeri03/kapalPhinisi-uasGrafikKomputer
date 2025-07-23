#include <iostream>
#include <vector>
#include <GL/freeglut.h>
#include <glm/glm.hpp>
#include <fstream> 

std::vector<glm::vec3> normals;

// Buffer Geometri Terpisah untuk Layar
std::vector<glm::vec3> sail_vertices;
std::vector<unsigned int> sail_indices;
std::vector<glm::vec3> sail_normals;
glm::vec3 bowsprit_end_point = {0.0f, 9.5f, 25.0f};
// Variabel Material Kayu (GLOBAL)
GLfloat wood_ambient[] = {0.4f, 0.25f, 0.15f, 1.0f};
GLfloat wood_diffuse[] = {0.6f, 0.4f, 0.2f, 1.0f};
GLfloat wood_specular[] = {0.5f, 0.5f, 0.5f, 1.0f};
GLfloat wood_shininess[] = {32.0f};

// Variabel Material Layar (GLOBAL)
GLfloat sail_ambient[] = {0.9f, 0.9f, 0.9f, 1.0f};
GLfloat sail_diffuse[] = {1.0f, 1.0f, 1.0f, 1.0f};
GLfloat sail_specular[] = {0.5f, 0.5f, 0.5f, 1.0f};
GLfloat sail_shininess[] = {50.0f};

// --- TITIK-TITIK IKAT UTAMA (sekarang global) ---
// Tiang Depan & Layanrya
glm::vec3 main_mast_pos = {0.0f, 7.5f, 5.0f};
float main_mast_height = 38.0f;
glm::vec3 main_p0, main_p1, main_p2, main_p3, main_top_p0, main_top_p1, main_top_p2;

// Tiang Belakang & Layarnya
glm::vec3 aft_mast_pos = {0.0f, 10.0f, -8.0f};
float aft_mast_height = 22.0f;
glm::vec3 aft_p0, aft_p1, aft_p2, aft_p3;


// =================================================================
// Variabel baru untuk Kontrol Kamera Mouse
float camera_distance = 40.0f;
float camera_angle_x = 25.0f;
float camera_angle_y = -45.0f;
glm::vec3 camera_target = {0.0f, 5.0f, 0.0f};
bool is_mouse_dragging = false;
int last_mouse_x = 0;
int last_mouse_y = 0;

void generateMastMesh(glm::vec3 position, float height, float radius, int sides);
void generateRailSegment(glm::vec3 p1, glm::vec3 p2, float radius);
// STRUKTUR DATA KURVA (TETAP SAMA)
struct HullCrossSection
{
    glm::vec2 P0, P1, P2, P3;
    glm::vec2 getPoint(float t) const {
        float u = 1.0f - t;
        float tt = t * t;
        float uu = u * u;
        float uuu = uu * u;
        float ttt = tt * t;
        return uuu * P0 + 3.0f * uu * t * P1 + 3.0f * u * tt * P2 + ttt * P3;
    }
};

// DEFINISI BLUEPRINT KURVA (TETAP SAMA)
HullCrossSection bowSection = { {0.0f, 0.0f}, {1.0f, 2.0f}, {3.0f, 7.5f}, {2.5f, 8.5f} };
HullCrossSection midSection = { {0.0f, 0.0f}, {4.5f, -1.5f}, {6.0f, 7.0f}, {5.5f, 8.0f} }; // P1.y diturunkan secara signifikan
HullCrossSection sternSection = { {0.0f, 0.0f}, {5.0f, -1.0f}, {7.0f, 6.0f}, {7.0f, 7.0f} };

// ===================================================================
// STRUKTUR DATA BARU UNTUK MESH 3D
// ===================================================================
std::vector<glm::vec3> vertices;
std::vector<unsigned int> indices;
float angle = 0.0f; // Untuk rotasi otomatis
const int segments = 20; // Jumlah titik per setengah kurva

// ===================================================================
// FUNGSI UTAMA UNTUK MEMBUAT MESH 3D
// ===================================================================
void generateHullMesh() {
    // 1. Hasilkan semua Vertices 3D
    // Z-posisi: buritan=-15, tengah=0, haluan=15
    std::vector<HullCrossSection> sections = {sternSection, midSection, bowSection};
    std::vector<float> z_positions = {-15.0f, 0.0f, 15.0f};
    int num_points_per_ring = 2 * segments;

    for (size_t i = 0; i < sections.size(); ++i) {
        // Hasilkan sisi kanan (x positif)
        for (int j = 0; j <= segments; ++j) {
            float t = static_cast<float>(j) / segments;
            glm::vec2 p2d = sections[i].getPoint(t);
            vertices.push_back(glm::vec3(p2d.x, p2d.y, z_positions[i]));
        }
        // Hasilkan sisi kiri (x negatif), lewati titik tengah (j=0)
        for (int j = segments - 1; j >= 1; --j) { // Perhatikan loop ini berjalan mundur
            float t = static_cast<float>(j) / segments;
            glm::vec2 p2d = sections[i].getPoint(t);
            vertices.push_back(glm::vec3(-p2d.x, p2d.y, z_positions[i]));
        }
    }

    // 2. Hasilkan Indices untuk "menjahit" kulit
    int num_points_per_section = 2 * segments; // Jumlah titik per irisan
    for (size_t i = 0; i < sections.size() - 1; ++i) {
        int section_start_index1 = i * num_points_per_section; // <-- PERBAIKAN 1
        int section_start_index2 = (i + 1) * num_points_per_section; // <-- PERBAIKAN 1

        for (int j = 0; j < num_points_per_section - 1; ++j) { // <-- PERBAIKAN 2
            // Bentuk segiempat (quad) dari 4 titik
            unsigned int p1 = section_start_index1 + j;
            unsigned int p2 = section_start_index1 + j + 1;
            unsigned int p3 = section_start_index2 + j;
            unsigned int p4 = section_start_index2 + j + 1;

            // Buat 2 segitiga dari segiempat tersebut
            indices.push_back(p1); indices.push_back(p3); indices.push_back(p4);
            indices.push_back(p1); indices.push_back(p4); indices.push_back(p2);
        }
        // Jahit quad terakhir untuk menutup celah
        int j = num_points_per_section - 1;
        unsigned int p1 = section_start_index1 + j;
        unsigned int p2 = section_start_index1 + 0; // Kembali ke titik awal
        unsigned int p3 = section_start_index2 + j;
        unsigned int p4 = section_start_index2 + 0; // Kembali ke titik awal

        indices.push_back(p1); indices.push_back(p3); indices.push_back(p4);
        indices.push_back(p1); indices.push_back(p4); indices.push_back(p2);
    }
    // Buat permukaan Dek dengan menghubungkan titik teratas
    for (size_t i = 0; i < sections.size() - 1; ++i) {
        int section_start_index1 = i * num_points_per_ring;
        int section_start_index2 = (i + 1) * num_points_per_ring;

        // Kita hanya akan menghubungkan separuh atas dari ring
        for (int j = segments / 2; j < num_points_per_ring - 1; ++j) {
            if (j > segments && j < segments + segments/2) continue; // Skip bagian bawah

            unsigned int p1 = section_start_index1 + j;
            unsigned int p2 = section_start_index1 + j + 1;
            unsigned int p3 = section_start_index2 + j;
            unsigned int p4 = section_start_index2 + j + 1;

            // Buat 2 segitiga, pastikan urutannya agar normal menghadap ke atas
            indices.push_back(p1); indices.push_back(p2); indices.push_back(p4);
            indices.push_back(p1); indices.push_back(p4); indices.push_back(p3);
        }
    }
    // ---------------------------------------------
   
    // Tutup bagian buritan (stern)
    int stern_start_index = 0; // Indeks awal dari vertex buritan
    

    for (int i = 1; i < num_points_per_ring - 1; ++i) {
        // Buat segitiga dari titik tengah (lunas), titik ke-i, dan titik ke-(i+1)
        indices.push_back(stern_start_index);
        indices.push_back(stern_start_index + i);
        indices.push_back(stern_start_index + i + 1);
    }
    // -------------------------------------------------
    // =================================================================
    // V TAMBAHKAN BLOK BARU INI UNTUK HALUAN MELENGKUNG V
    // =================================================================
    // Buat haluan melengkung dengan menyatukan semua titik ke satu titik depan
    int bow_start_index = 2 * num_points_per_ring;
        
    // Tambahkan satu vertex baru sebagai ujung paling depan (stem)
    glm::vec3 stem_point = bowsprit_end_point; // Posisinya di depan & tengah
    int stem_index = vertices.size();
    vertices.push_back(stem_point);
        
    // Hubungkan semua titik di cincin haluan ke titik stem
    for (int i = 0; i < num_points_per_ring; ++i) {
        indices.push_back(bow_start_index + i);
        // Gunakan modulo untuk menghubungkan titik terakhir kembali ke awal
        indices.push_back(bow_start_index + ((i + 1) % num_points_per_ring));
        indices.push_back(stem_index);
    }
    // =================================================================
}

// GANTI SELURUH FUNGSI DECKHOUSE ANDA DENGAN INI
void generateDeckhouseMesh() {
    // Helper lambda untuk membuat kotak dengan mudah
    auto generateBox = [&](glm::vec3 pos, glm::vec3 size) {
        int base_index = vertices.size();
        vertices.push_back(pos); // 0
        vertices.push_back(pos + glm::vec3(size.x, 0, 0)); // 1
        vertices.push_back(pos + glm::vec3(size.x, size.y, 0)); // 2
        vertices.push_back(pos + glm::vec3(0, size.y, 0)); // 3
        vertices.push_back(pos + glm::vec3(0, 0, size.z)); // 4
        vertices.push_back(pos + glm::vec3(size.x, 0, size.z)); // 5
        vertices.push_back(pos + glm::vec3(size.x, size.y, size.z)); // 6
        vertices.push_back(pos + glm::vec3(0, size.y, size.z)); // 7
        
        int indices_data[] = {0,2,1, 0,3,2, 1,6,5, 1,2,6, 4,3,0, 4,7,3, 5,7,4, 5,6,7, 3,6,7, 3,2,6, 4,1,5, 4,0,1};
        for(int i=0; i<36; ++i) indices.push_back(base_index + indices_data[i]);
    };

    // --- 1. Kabin Utama Tingkat Bawah ---
    glm::vec3 lower_cabin_pos = {-4.5f, 7.5f, -14.0f};
    glm::vec3 lower_cabin_size = {9.0f, 3.0f, 18.0f};
    generateBox(lower_cabin_pos, lower_cabin_size);

    // --- 2. Lantai Dek Atas ---
    glm::vec3 upper_deck_pos = {-5.0f, 10.5f, -14.5f};
    glm::vec3 upper_deck_size = {10.0f, 0.2f, 19.0f};
    generateBox(upper_deck_pos, upper_deck_size);

    // --- 3. Kabin Tingkat Atas (Jembatan Anjungan) ---
    glm::vec3 upper_cabin_pos = {-3.5f, 10.7f, -12.0f};
    glm::vec3 upper_cabin_size = {7.0f, 2.5f, 8.0f};
    generateBox(upper_cabin_pos, upper_cabin_size);

    // --- 4. Jendela-jendela Kabin Bawah (Sebagai kotak hitam kecil) ---
    float window_size = 0.8f;
    float window_depth = 0.1f;
    for (float z = -12.0f; z < 3.0f; z += 2.5f) {
        // Jendela kanan
        generateBox({lower_cabin_pos.x + lower_cabin_size.x, 8.5f, z}, {window_depth, window_size, window_size});
        // Jendela kiri
        generateBox({lower_cabin_pos.x - window_depth, 8.5f, z}, {window_depth, window_size, window_size});
    }
    
    // --- 5. Pagar di Dek Atas ---
    float rail_height = 1.2f;
    float rail_radius = 0.05f;
    float rail_post_radius = 0.08f;
    // Pagar sisi kanan
    for(float z = -14.0f; z < 4.0f; z += 2.0f) {
        generateMastMesh({upper_deck_pos.x + upper_deck_size.x, upper_deck_pos.y, z}, rail_height, rail_post_radius, 8);
    }
    generateRailSegment({upper_deck_pos.x + upper_deck_size.x, upper_deck_pos.y + rail_height, -14.0f},
                          {upper_deck_pos.x + upper_deck_size.x, upper_deck_pos.y + rail_height, 4.0f}, rail_radius);
    // Pagar sisi kiri
    for(float z = -14.0f; z < 4.0f; z += 2.0f) {
        generateMastMesh({upper_deck_pos.x, upper_deck_pos.y, z}, rail_height, rail_post_radius, 8);
    }
    generateRailSegment({upper_deck_pos.x, upper_deck_pos.y + rail_height, -14.0f},
                          {upper_deck_pos.x, upper_deck_pos.y + rail_height, 4.0f}, rail_radius);
    // Pagar belakang
    generateRailSegment({upper_deck_pos.x, upper_deck_pos.y + rail_height, -14.0f},
                          {upper_deck_pos.x + upper_deck_size.x, upper_deck_pos.y + rail_height, -14.0f}, rail_radius);
}

void generateMastMesh(glm::vec3 position, float height, float radius, int sides) {
    #define M_PI 3.14159265358979323846 // Pastikan M_PI terdefinisi

    int base_index = vertices.size();

    // Buat vertex untuk sisi silinder
    for (int i = 0; i < sides; ++i) {
        float angle = 2.0f * M_PI * static_cast<float>(i) / sides;
        float x = cos(angle) * radius;
        float z = sin(angle) * radius;

        // Tambahkan vertex bawah dan atas untuk setiap segmen
        vertices.push_back(position + glm::vec3(x, 0, z));      // Vertex bawah
        vertices.push_back(position + glm::vec3(x, height, z)); // Vertex atas
    }

    // Buat indices untuk sisi silinder
for (int i = 0; i < sides; ++i) {
    // Tentukan 4 titik untuk membuat segiempat
    int p1_bottom = base_index + i * 2;
    int p1_top = base_index + i * 2 + 1;
    // Gunakan modulo (%) untuk menyambung kembali ke titik awal
    int p2_bottom = base_index + ((i + 1) % sides) * 2;
    int p2_top = base_index + ((i + 1) % sides) * 2 + 1;

    // Segitiga pertama
    indices.push_back(p1_bottom); indices.push_back(p2_bottom); indices.push_back(p2_top);
    // Segitiga kedua
    indices.push_back(p1_bottom); indices.push_back(p2_top); indices.push_back(p1_top);
}

    // Buat tutup atas dan bawah (Triangle Fan)
    // Pertama, buat vertex pusat untuk tutup atas dan bawah
    vertices.push_back(position + glm::vec3(0, height, 0)); // Pusat atas
    int top_center_index = vertices.size() - 1;
    vertices.push_back(position + glm::vec3(0, 0, 0)); // Pusat bawah
    int bottom_center_index = vertices.size() - 1;

    for (int i = 0; i < sides; ++i) {
        int p1_top = base_index + i * 2 + 1;
        int p2_top = base_index + ((i + 1) % sides) * 2 + 1;
        // Tutup atas
        indices.push_back(top_center_index); indices.push_back(p2_top); indices.push_back(p1_top);

        int p1_bottom = base_index + i * 2;
        int p2_bottom = base_index + ((i + 1) % sides) * 2;
        // Tutup bawah
        indices.push_back(bottom_center_index); indices.push_back(p1_bottom); indices.push_back(p2_bottom);
    }
}

void recalculateAllNormals() {
    normals.assign(vertices.size(), glm::vec3(0.0f, 0.0f, 0.0f));
    // Hitung normal untuk setiap segitiga dan akumulasikan
    for (size_t i = 0; i < indices.size(); i += 3) {
        glm::vec3 p1 = vertices[indices[i]];
        glm::vec3 p2 = vertices[indices[i+1]];
        glm::vec3 p3 = vertices[indices[i+2]];

        glm::vec3 v1 = p2 - p1;
        glm::vec3 v2 = p3 - p1;
        glm::vec3 normal = glm::cross(v1, v2);

        normals[indices[i]] += normal;
        normals[indices[i+1]] += normal;
        normals[indices[i+2]] += normal;
    }

    // Normalisasi semua normal (buat panjangnya menjadi 1)
    for (size_t i = 0; i < normals.size(); ++i) {
        normals[i] = glm::normalize(normals[i]);
    }
}

// FUNGSI HELPER BARU UNTUK MENDAPATKAN LEBAR DEK SECARA DINAMIS
float getDeckEdgeXAtZ(float z) {
    // Tentukan di antara seksi mana posisi Z berada
    if (z >= -15.0f && z < 0.0f) {
        // Antara Buritan (stern) dan Tengah (mid)
        float t = (z - (-15.0f)) / (0.0f - (-15.0f)); // Faktor interpolasi 0..1
        // Interpolasi antara P3.x buritan (7.0) dan P3.x tengah (5.5)
        return glm::mix(sternSection.P3.x, midSection.P3.x, t);
    } else if (z >= 0.0f && z <= 15.0f) {
        // Antara Tengah (mid) dan Haluan (bow)
        float t = (z - 0.0f) / (15.0f - 0.0f); // Faktor interpolasi 0..1
        // Interpolasi antara P3.x tengah (5.5) dan P3.x haluan (2.5)
        return glm::mix(midSection.P3.x, bowSection.P3.x, t);
    }
    return 5.0f; // Default fallback
}

// FUNGSI HELPER BARU UNTUK MENDAPATKAN KETINGGIAN DEK SECARA DINAMIS
float getDeckEdgeYAtZ(float z) {
    if (z >= -15.0f && z < 0.0f) {
        // Antara Buritan (stern) dan Tengah (mid)
        float t = (z - (-15.0f)) / (0.0f - (-15.0f));
        // Interpolasi antara P3.y buritan (7.0) dan P3.y tengah (8.0)
        return glm::mix(sternSection.P3.y, midSection.P3.y, t);
    } else if (z >= 0.0f && z <= 15.0f) {
        // Antara Tengah (mid) dan Haluan (bow)
        float t = (z - 0.0f) / (15.0f - 0.0f);
        // Interpolasi antara P3.y tengah (8.0) dan P3.y haluan (8.5)
        return glm::mix(midSection.P3.y, bowSection.P3.y, t);
    }
    return 7.5f; // Default fallback
}

// FUNGSI BARU UNTUK MEMBUAT SEGMEN REL
void generateRailSegment(glm::vec3 p1, glm::vec3 p2, float radius) {
    int base_index = vertices.size();

    glm::vec3 up(0.0f, 1.0f, 0.0f);
    glm::vec3 dir = glm::normalize(p2 - p1);
    glm::vec3 side = glm::normalize(glm::cross(dir, up)) * radius;
    up = glm::normalize(glm::cross(side, dir)) * radius;

    // 8 titik sudut dari balok rel
    vertices.push_back(p1 - side - up);
    vertices.push_back(p1 + side - up);
    vertices.push_back(p1 + side + up);
    vertices.push_back(p1 - side + up);
    vertices.push_back(p2 - side - up);
    vertices.push_back(p2 + side - up);
    vertices.push_back(p2 + side + up);
    vertices.push_back(p2 - side + up);

    // Indices untuk 6 sisi balok
    int indices_data[] = {
        0, 1, 2, 0, 2, 3, // Depan
        4, 5, 6, 4, 6, 7, // Belakang
        3, 2, 6, 3, 6, 7, // Atas
        0, 1, 5, 0, 5, 4, // Bawah
        4, 7, 3, 4, 3, 0, // Kiri
        1, 5, 6, 1, 6, 2  // Kanan
    };
    
    for(int i = 0; i < 36; ++i) {
        indices.push_back(base_index + indices_data[i]);
    }
}
void generateRailings() {
    float post_height = 1.5f;
    float post_radius = 0.08f;
    int post_sides = 8;
    float rail_radius = 0.05f;

    // --- Tahap 1: Tentukan semua titik jalur pagar ---
    std::vector<glm::vec3> right_path;
    std::vector<glm::vec3> left_path;

    // Buat jalur untuk sisi kanan dan kiri
    for (float z = -15.0f; z <= 15.0f; z += 2.0f) {
        float x_pos = getDeckEdgeXAtZ(z);
        float y_pos = getDeckEdgeYAtZ(z);
        right_path.push_back({x_pos, y_pos, z});
        left_path.push_back({-x_pos, y_pos, z});
    }

    // --- Tahap 2: Buat mesh berdasarkan jalur yang sudah ada ---
    
    // Buat tiang dan rel sisi kanan
    for (size_t i = 0; i < right_path.size(); ++i) {
        generateMastMesh(right_path[i], post_height, post_radius, post_sides);
        if (i < right_path.size() - 1) {
            generateRailSegment(right_path[i] + glm::vec3(0, post_height, 0),
                                  right_path[i+1] + glm::vec3(0, post_height, 0), rail_radius);
        }
    }

    // Buat tiang dan rel sisi kiri
    for (size_t i = 0; i < left_path.size(); ++i) {
        generateMastMesh(left_path[i], post_height, post_radius, post_sides);
        if (i < left_path.size() - 1) {
            generateRailSegment(left_path[i] + glm::vec3(0, post_height, 0),
                                  left_path[i+1] + glm::vec3(0, post_height, 0), rail_radius);
        }
    }

    // Buat tiang dan rel pagar belakang
    float back_y = getDeckEdgeYAtZ(-15.0f);
    float back_x_extent = getDeckEdgeXAtZ(-15.0f);
    for (float x = -back_x_extent; x <= back_x_extent; x += 1.5f) {
        generateMastMesh({x, back_y, -15.0f}, post_height, post_radius, post_sides);
    }
    generateRailSegment({-back_x_extent, back_y + post_height, -15.0f},
                          {back_x_extent, back_y + post_height, -15.0f}, rail_radius);
    
    float front_z = 15.0f;
    float front_y = getDeckEdgeYAtZ(front_z);
    float front_x_extent = getDeckEdgeXAtZ(front_z);
    for (float x = -front_x_extent; x <= front_x_extent; x += 1.0f) { // Jarak antar tiang lebih rapat di depan
        generateMastMesh({x, front_y, front_z}, post_height, post_radius, post_sides);
    }
    // Hubungkan rel depan dengan rel samping
    // Titik terakhir dari jalur samping adalah titik depan
    glm::vec3 front_right_top = right_path.back() + glm::vec3(0, post_height, 0);
    glm::vec3 front_left_top = left_path.back() + glm::vec3(0, post_height, 0);
    generateRailSegment(front_left_top, front_right_top, rail_radius);
}

// FUNGSI BARU UNTUK MEMBUAT LAYAR
void generateSailMesh(glm::vec3 p0, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3, int divisions, 
                      std::vector<glm::vec3>& target_vertices, 
                      std::vector<unsigned int>& target_indices) {
    int base_index = vertices.size();
    for (int i = 0; i <= divisions; ++i) {
        float t = static_cast<float>(i) / divisions;
        glm::vec3 side1 = glm::mix(p0, p3, t);
        glm::vec3 side2 = glm::mix(p1, p2, t);
        for (int j = 0; j <= divisions; ++j) {
            float s = static_cast<float>(j) / divisions;
            glm::vec3 point = glm::mix(side1, side2, s);
            // Deformasi dinonaktifkan untuk tampilan 2D
            // Deformasi untuk efek "menggembung" (DIAKTIFKAN KEMBALI)
            // glm::vec3 sail_normal = glm::normalize(glm::cross(p1 - p0, p3 - p0));
            // float billow = sin(t * M_PI) * sin(s * M_PI) * 2.5f; // Kita kurangi sedikit lengkungannya
            // point += sail_normal * billow;
            vertices.push_back(point);

        }
    }
    for (int i = 0; i < divisions; ++i) {
        for (int j = 0; j < divisions; ++j) {
            int row1 = i * (divisions + 1);
            int row2 = (i + 1) * (divisions + 1);
            indices.push_back(base_index + row1 + j);
            indices.push_back(base_index + row1 + j + 1);
            indices.push_back(base_index + row2 + j + 1);
            indices.push_back(base_index + row1 + j);
            indices.push_back(base_index + row2 + j + 1);
            indices.push_back(base_index + row2 + j);
        }
    }
}

void generateBowsprit() {
    // --- 1. Definisikan Titik & Parameter Utama ---
    glm::vec3 base_point = {0.0f, 8.5f, 15.0f}; // Pangkal di dek, tetap
    // Ujung bowsprit sekarang menggunakan variabel global
    glm::vec3 end_point = bowsprit_end_point; 
    
    float main_spar_radius = 0.2f;
    float support_post_height = 1.2f;
    float support_post_radius = 0.07f;
    float rail_radius = 0.05f;

    // --- 2. Buat Dua Tiang Penyangga Utama (Bentuk V) ---
    // Ini akan menjadi fondasi struktur baru kita
    generateRailSegment(base_point + glm::vec3(-0.5f, 0, 0), end_point, main_spar_radius);
    generateRailSegment(base_point + glm::vec3(0.5f, 0, 0), end_point, main_spar_radius);

    // --- 3. Buat Tiang-tiang Vertikal (Stanchions) di Atas Bowsprit ---
    std::vector<glm::vec3> stanchion_tops;
    for (int i = 1; i <= 3; ++i) {
        float t = i / 3.5f; // Posisi di sepanjang bowsprit
        glm::vec3 base_on_spar = glm::mix(base_point, end_point, t);

        // Tiang kanan
        glm::vec3 post_base_right = base_on_spar + glm::vec3(0.5f + t*0.2f, main_spar_radius, 0);
        generateMastMesh(post_base_right, support_post_height, support_post_radius, 8);
        stanchion_tops.push_back(post_base_right + glm::vec3(0, support_post_height, 0));

        // Tiang kiri
        glm::vec3 post_base_left = base_on_spar + glm::vec3(-0.5f - t*0.2f, main_spar_radius, 0);
        generateMastMesh(post_base_left, support_post_height, support_post_radius, 8);
        stanchion_tops.push_back(post_base_left + glm::vec3(0, support_post_height, 0));
    }

    // --- 4. Buat Pagar Depan (Pulpit) dengan Menghubungkan Tiang Vertikal ---
    // Hubungkan sisi kanan
    generateRailSegment(stanchion_tops[0], stanchion_tops[2], rail_radius);
    generateRailSegment(stanchion_tops[2], stanchion_tops[4], rail_radius);
    // Hubungkan sisi kiri
    generateRailSegment(stanchion_tops[1], stanchion_tops[3], rail_radius);
    generateRailSegment(stanchion_tops[3], stanchion_tops[5], rail_radius);
    // Hubungkan bagian paling depan
    generateRailSegment(stanchion_tops[4], stanchion_tops[5], rail_radius);
}

// FUNGSI BARU KHUSUS UNTUK TALI-TEMALI KOMPLEKS
// GANTI FUNGSI LAMA DENGAN VERSI FINAL INI
void generateRigging() {
    float rope_radius = 0.04f;
    float thin_rope_radius = 0.02f;

    // Ambil Titik-titik Kunci dari Variabel Global
    glm::vec3 main_mast_top = main_mast_pos + glm::vec3(0, main_mast_height, 0);
    glm::vec3 aft_mast_top = aft_mast_pos + glm::vec3(0, aft_mast_height, 0);

    // ======================================================
    // --- RIGGING TIANG DEPAN (MAIN MAST) ---
    // ======================================================
    
    // 1. Forestay Utama (ke ujung bowsprit)
    generateRailSegment(main_mast_top, bowsprit_end_point, rope_radius);

    // 2. Shrouds (Tali Samping)
    for (int i = 0; i < 4; ++i) {
        float z_pos = -2.0f + (i * 4.0f);
        float x_pos = getDeckEdgeXAtZ(z_pos) + 0.5f;
        float y_pos = getDeckEdgeYAtZ(z_pos);
        generateRailSegment(main_mast_top, {x_pos, y_pos, z_pos}, rope_radius);
        generateRailSegment(main_mast_top, {-x_pos, y_pos, z_pos}, rope_radius);
    }

    // DIHAPUS: Backstay ke tiang belakang yang membuat berantakan

    // ======================================================
    // --- RIGGING TIANG BELAKANG (AFT MAST) ---
    // ======================================================

    // 1. Shrouds (Tali Samping)
    for (int i = 0; i < 3; ++i) {
        float z_pos = -13.0f + (i * 3.0f);
        float x_pos = getDeckEdgeXAtZ(z_pos) + 0.5f;
        float y_pos = getDeckEdgeYAtZ(z_pos);
        generateRailSegment(aft_mast_top, {x_pos, y_pos, z_pos}, rope_radius);
        generateRailSegment(aft_mast_top, {-x_pos, y_pos, z_pos}, rope_radius);
    }
    
    // 2. Backstay (ke buritan)
    generateRailSegment(aft_mast_top, {0.0f, getDeckEdgeYAtZ(-15.0f), -15.0f}, rope_radius);

    // ======================================================
    // --- TALI LAYAR (SHEETS) ---
    // ======================================================
    
    // 1. Tali Pengikat Ujung Layar (Sheets)
    // Untuk layar utama depan - Diikat ke belakang kabin
    generateRailSegment(main_p1, {getDeckEdgeXAtZ(-4.0f), getDeckEdgeYAtZ(-4.0f) + 1.5f, -4.0f}, thin_rope_radius);
    // Untuk layar utama belakang - Diikat ke pagar belakang
    generateRailSegment(aft_p1, {getDeckEdgeXAtZ(-15.0f), getDeckEdgeYAtZ(-15.0f) + 1.5f, -15.0f}, thin_rope_radius);

    // DIHAPUS: Topping lifts yang membuat berantakan
}

void generateLadders() {
    // --- Parameter Tangga ---
    float ladder_width = 0.6f;
    float rail_radius = 0.05f;  // Ketebalan pegangan samping
    float rung_radius = 0.03f;  // Ketebalan anak tangga

    // ======================================================
    // --- TANGGA TIANG DEPAN (DIPINDAH KE SISI KANAN) ---
    // ======================================================
    float main_mast_radius = 0.5f;
    float main_ladder_bottom_y = main_mast_pos.y + 2.0f;
    float main_ladder_top_y = main_mast_pos.y + main_mast_height * 0.9f;

    // --- Bagian 1: Tangga di Sisi Kanan (Sudah Ada) ---
    float side_offset_x_right = main_mast_radius + rail_radius;
    // Buat pegangan vertikal
    generateRailSegment({main_mast_pos.x + side_offset_x_right, main_ladder_bottom_y, main_mast_pos.z - ladder_width/2}, {main_mast_pos.x + side_offset_x_right, main_ladder_top_y, main_mast_pos.z - ladder_width/2}, rail_radius);
    generateRailSegment({main_mast_pos.x + side_offset_x_right, main_ladder_bottom_y, main_mast_pos.z + ladder_width/2}, {main_mast_pos.x + side_offset_x_right, main_ladder_top_y, main_mast_pos.z + ladder_width/2}, rail_radius);
    // Buat anak tangga
    for (float y = main_ladder_bottom_y; y < main_ladder_top_y; y += 0.4f) {
        glm::vec3 start = {main_mast_pos.x + side_offset_x_right, y, main_mast_pos.z - ladder_width/2};
        glm::vec3 end   = {main_mast_pos.x + side_offset_x_right, y, main_mast_pos.z + ladder_width/2};
        generateRailSegment(start, end, rung_radius);
    }

    // --- Bagian 2: Tangga Duplikat di Sisi Kiri (Baru) ---
    float side_offset_x_left = -main_mast_radius - rail_radius;
    // Buat pegangan vertikal
    generateRailSegment({main_mast_pos.x + side_offset_x_left, main_ladder_bottom_y, main_mast_pos.z - ladder_width/2}, {main_mast_pos.x + side_offset_x_left, main_ladder_top_y, main_mast_pos.z - ladder_width/2}, rail_radius);
    generateRailSegment({main_mast_pos.x + side_offset_x_left, main_ladder_bottom_y, main_mast_pos.z + ladder_width/2}, {main_mast_pos.x + side_offset_x_left, main_ladder_top_y, main_mast_pos.z + ladder_width/2}, rail_radius);
    // Buat anak tangga
    for (float y = main_ladder_bottom_y; y < main_ladder_top_y; y += 0.4f) {
        glm::vec3 start = {main_mast_pos.x + side_offset_x_left, y, main_mast_pos.z - ladder_width/2};
        glm::vec3 end   = {main_mast_pos.x + side_offset_x_left, y, main_mast_pos.z + ladder_width/2};
        generateRailSegment(start, end, rung_radius);
    }

    // ======================================================
    // --- TANGGA TIANG BELAKANG (DI DEPAN TIANG) ---
    // ======================================================
    float aft_mast_radius = 0.4f;
    float aft_ladder_bottom_y = aft_mast_pos.y + 2.0f;
    float aft_ladder_top_y = aft_mast_pos.y + aft_mast_height * 0.9f;

    // 1. Tentukan posisi untuk dua pegangan vertikal di depan tiang
    glm::vec3 aft_left_rail_bottom = {aft_mast_pos.x - ladder_width/2, aft_ladder_bottom_y, aft_mast_pos.z + aft_mast_radius};
    glm::vec3 aft_left_rail_top    = {aft_mast_pos.x - ladder_width/2, aft_ladder_top_y,    aft_mast_pos.z + aft_mast_radius};

    glm::vec3 aft_right_rail_bottom = {aft_mast_pos.x + ladder_width/2, aft_ladder_bottom_y, aft_mast_pos.z + aft_mast_radius};
    glm::vec3 aft_right_rail_top    = {aft_mast_pos.x + ladder_width/2, aft_ladder_top_y,    aft_mast_pos.z + aft_mast_radius};

    // 2. Buat dua pegangan vertikal
    generateRailSegment(aft_left_rail_bottom, aft_left_rail_top, rail_radius);
    generateRailSegment(aft_right_rail_bottom, aft_right_rail_top, rail_radius);

    // 3. Buat anak tangga horizontal
    for (float y = aft_ladder_bottom_y; y < aft_ladder_top_y; y += 0.4f) {
        glm::vec3 start = {aft_mast_pos.x - ladder_width/2, y, aft_mast_pos.z + aft_mast_radius};
        glm::vec3 end   = {aft_mast_pos.x + ladder_width/2, y, aft_mast_pos.z + aft_mast_radius};
        generateRailSegment(start, end, rung_radius);
    }
}

// --- FUNGSI BARU UNTUK MEMBUAT TANGGA KABIN ---
void generateCabinLadder() {
    // --- Parameter Tangga ---
    float ladder_width = 1.0f;      // Lebar anak tangga
    int num_steps = 8;              // Jumlah anak tangga
    float rail_radius = 0.05f;      // Ketebalan pegangan
    float rung_thickness = 0.05f;   // Ketebalan papan anak tangga

    // --- Posisi Tangga ---
    // Posisi ini diatur agar menempel di sisi kiri kabin bawah
    // dan mengarah miring ke bawah dari dek atas ke dek utama.
    glm::vec3 bottom_pos = {-2.0f, 7.7f, -13.0f}; // Titik awal di dek utama belakang
    glm::vec3 top_pos = {-2.0f, 10.7f, -12.0f};// Titik akhir di dek bawah

    // --- 1. Buat Pegangan (Railings) Samping ---
    // Pegangan kiri
    generateRailSegment(top_pos, bottom_pos, rail_radius);
    // Pegangan kanan (geser sejauh lebar tangga)
    generateRailSegment(top_pos + glm::vec3(0, 0, ladder_width), 
                        bottom_pos + glm::vec3(0, 0, ladder_width), rail_radius);

    // --- 2. Buat Anak Tangga ---
    // Gunakan loop untuk menempatkan setiap anak tangga secara merata
    // dari posisi atas ke bawah.
    for (int i = 0; i <= num_steps; ++i) {
        // Hitung faktor interpolasi (t) dari 0.0 sampai 1.0
        float t = static_cast<float>(i) / num_steps;

        // Dapatkan posisi tengah anak tangga dengan interpolasi linear
        glm::vec3 step_center = glm::mix(top_pos, bottom_pos, t);

        // Definisikan titik awal dan ukuran untuk anak tangga (sebagai balok)
        glm::vec3 step_start_pos = step_center - glm::vec3(0, rung_thickness / 2.0f, 0);
        glm::vec3 step_size = {0.0f, rung_thickness, ladder_width}; // Balok tipis

        // Gunakan fungsi helper 'generateBox' yang sudah ada
        int base_index = vertices.size();
        vertices.push_back(step_start_pos); // 0
        vertices.push_back(step_start_pos + glm::vec3(step_size.x, 0, 0)); // 1
        vertices.push_back(step_start_pos + glm::vec3(step_size.x, step_size.y, 0)); // 2
        vertices.push_back(step_start_pos + glm::vec3(0, step_size.y, 0)); // 3
        vertices.push_back(step_start_pos + glm::vec3(0, 0, step_size.z)); // 4
        vertices.push_back(step_start_pos + glm::vec3(step_size.x, 0, step_size.z)); // 5
        vertices.push_back(step_start_pos + glm::vec3(step_size.x, step_size.y, step_size.z)); // 6
        vertices.push_back(step_start_pos + glm::vec3(0, step_size.y, step_size.z)); // 7
        
        int indices_data[] = {0,2,1, 0,3,2, 1,6,5, 1,2,6, 4,3,0, 4,7,3, 5,7,4, 5,6,7, 3,6,7, 3,2,6, 4,1,5, 4,0,1};
        for(int j=0; j<36; ++j) indices.push_back(base_index + indices_data[j]);
    }
}

void recalculateSailNormals() {
    sail_normals.assign(sail_vertices.size(), glm::vec3(0.0f, 0.0f, 0.0f));
    for (size_t i = 0; i < sail_indices.size(); i += 3) {
        glm::vec3 p1 = sail_vertices[sail_indices[i]];
        glm::vec3 p2 = sail_vertices[sail_indices[i+1]];
        glm::vec3 p3 = sail_vertices[sail_indices[i+2]];
        glm::vec3 v1 = p2 - p1;
        glm::vec3 v2 = p3 - p1;
        glm::vec3 normal = glm::cross(v1, v2);
        sail_normals[sail_indices[i]] += normal;
        sail_normals[sail_indices[i+1]] += normal;
        sail_normals[sail_indices[i+2]] += normal;
    }
    for (size_t i = 0; i < sail_normals.size(); ++i) {
        sail_normals[i] = glm::normalize(sail_normals[i]);
    }
}
//=========KONTROL KAMERA MENGGUNAKAN MOUSE==============
// Fungsi ini dipanggil saat tombol mouse ditekan atau dilepas
void mouseButton(int button, int state, int x, int y) {
    if (button == GLUT_LEFT_BUTTON) {
        if (state == GLUT_DOWN) {
            is_mouse_dragging = true;
            last_mouse_x = x;
            last_mouse_y = y;
        } else {
            is_mouse_dragging = false;
        }
    }
}

// Fungsi ini dipanggil saat mouse bergerak sambil menahan tombol
void mouseMove(int x, int y) {
    if (is_mouse_dragging) {
        float dx = x - last_mouse_x;
        float dy = y - last_mouse_y;

        camera_angle_y += dx * 0.5f; // Kontrol rotasi horizontal
        camera_angle_x += dy * 0.5f; // Kontrol rotasi vertikal

        // Batasi sudut vertikal agar tidak terbalik
        if (camera_angle_x > 89.0f) camera_angle_x = 89.0f;
        if (camera_angle_x < -89.0f) camera_angle_x = -89.0f;

        last_mouse_x = x;
        last_mouse_y = y;

        glutPostRedisplay(); // Minta untuk menggambar ulang
    }
}

// Fungsi ini dipanggil saat scroll wheel digunakan (untuk zoom)
void mouseWheel(int button, int dir, int x, int y) {
    if (dir > 0) {
        camera_distance -= 2.0f; // Zoom in
    } else {
        camera_distance += 2.0f; // Zoom out
    }

    // Batasi jarak zoom
    if (camera_distance < 5.0f) camera_distance = 5.0f;
    if (camera_distance > 100.0f) camera_distance = 100.0f;

    glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y) {
    float pan_speed = 0.5f; // Kecepatan gerak kamera

    switch (key) {
        case 'w': // Geser maju (sepanjang sumbu Z)
            camera_target.z -= pan_speed;
            break;
        case 's': // Geser mundur (sepanjang sumbu Z)
            camera_target.z += pan_speed;
            break;
        case 'a': // Geser kiri (sepanjang sumbu X)
            camera_target.x -= pan_speed;
            break;
        case 'd': // Geser kanan (sepanjang sumbu X)
            camera_target.x += pan_speed;
            break;
        case 'q': // Geser atas (sepanjang sumbu Y)
            camera_target.y += pan_speed;
            break;
        case 'e': // Geser bawah (sepanjang sumbu Y)
            camera_target.y -= pan_speed;
            break;
    }
    glutPostRedisplay(); // Minta untuk menggambar ulang scene
}

// --- FUNGSI BARU UNTUK EKSPOR KE FILE .OBJ ---
void exportToObj(const char* filename) {
    // Buka file untuk ditulis
    std::ofstream objFile(filename);
    if (!objFile.is_open()) {
        std::cerr << "Error: Tidak bisa membuka file untuk ekspor!" << std::endl;
        return;
    }

    objFile << "# Model Kapal Phinisi diekspor dari OpenGL" << std::endl;
    objFile << "# Vertices: " << vertices.size() << std::endl;
    objFile << "# Normals: " << normals.size() << std::endl;
    objFile << "# Faces: " << indices.size() / 3 << std::endl;

    // 1. Tulis semua data vertex (v)
    for (const auto& vertex : vertices) {
        objFile << "v " << vertex.x << " " << vertex.y << " " << vertex.z << std::endl;
    }

    // 2. Tulis semua data normal (vn)
    for (const auto& normal : normals) {
        objFile << "vn " << normal.x << " " << normal.y << " " << normal.z << std::endl;
    }

    // 3. Tulis semua data face/indices (f)
    // Format .obj adalah 1-based, jadi kita perlu menambahkan 1 pada setiap index
    for (size_t i = 0; i < indices.size(); i += 3) {
        unsigned int i1 = indices[i] + 1;
        unsigned int i2 = indices[i+1] + 1;
        unsigned int i3 = indices[i+2] + 1;
        // Format: f v1//vn1 v2//vn2 v3//vn3 (tanpa texture coordinates)
        objFile << "f " << i1 << "//" << i1 << " " << i2 << "//" << i2 << " " << i3 << "//" << i3 << std::endl;
    }
    
    // Lakukan hal yang sama untuk layar jika ada
    if (!sail_vertices.empty()) {
        int vertex_offset = vertices.size();
        for (const auto& vertex : sail_vertices) {
            objFile << "v " << vertex.x << " " << vertex.y << " " << vertex.z << std::endl;
        }
        for (const auto& normal : sail_normals) {
            objFile << "vn " << normal.x << " " << normal.y << " " << normal.z << std::endl;
        }
        for (size_t i = 0; i < sail_indices.size(); i += 3) {
            unsigned int i1 = vertex_offset + sail_indices[i] + 1;
            unsigned int i2 = vertex_offset + sail_indices[i+1] + 1;
            unsigned int i3 = vertex_offset + sail_indices[i+2] + 1;
            objFile << "f " << i1 << "//" << i1 << " " << i2 << "//" << i2 << " " << i3 << "//" << i3 << std::endl;
        }
    }


    objFile.close();
    std::cout << "Model berhasil diekspor ke " << filename << std::endl;
}
// ===================================================================
// FUNGSI DISPLAY OPENGL (KINI DALAM 3D)
// ===================================================================
void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0, 1.0, 1.0, 100.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // Hitung posisi kamera berdasarkan sudut dan jarak DARI TARGET
    float rad_x = camera_angle_x * M_PI / 180.0f;
    float rad_y = camera_angle_y * M_PI / 180.0f;

    // Hitung offset posisi mata dari target
    float eye_offset_x = camera_distance * cos(rad_x) * sin(rad_y);
    float eye_offset_y = camera_distance * sin(rad_x);
    float eye_offset_z = camera_distance * cos(rad_x) * cos(rad_y);

    // Posisi mata adalah posisi target ditambah offset
    glm::vec3 eye_pos = camera_target + glm::vec3(eye_offset_x, eye_offset_y, eye_offset_z);

    // Atur kamera untuk melihat dari eye_pos ke camera_target
    gluLookAt(eye_pos.x, eye_pos.y, eye_pos.z,
              camera_target.x, camera_target.y, camera_target.z,
              0.0, 1.0, 0.0);
    
    // Atur properti pencahayaan
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_NORMALIZE);
    GLfloat light_pos[] = {-25.0f, 30.0f, 30.0f, 1.0f};
    glLightfv(GL_LIGHT0, GL_POSITION, light_pos);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    // --- GAMBAR LAMBUNG & STRUKTUR KAYU ---
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, wood_ambient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, wood_diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, wood_specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, wood_shininess);
    
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glVertexPointer(3, GL_FLOAT, 0, vertices.data());
    glNormalPointer(GL_FLOAT, 0, normals.data());
    glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, indices.data());
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_VERTEX_ARRAY);

    // --- GAMBAR LAYAR ---
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, sail_ambient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, sail_diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, sail_specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, sail_shininess);

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glVertexPointer(3, GL_FLOAT, 0, sail_vertices.data());
    glNormalPointer(GL_FLOAT, 0, sail_normals.data());
    glDrawElements(GL_TRIANGLES, sail_indices.size(), GL_UNSIGNED_INT, sail_indices.data());
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_VERTEX_ARRAY);

    glutSwapBuffers();
}

// Timer untuk animasi rotasi
void timer(int) {
    angle += 0.5f; // Kecepatan rotasi
    if (angle > 360.0f) angle -= 360.0f;
    glutPostRedisplay();
    glutTimerFunc(16, timer, 0); // ~60 FPS
}

// FUNGSI MAIN
int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(800, 800);
    glutCreateWindow("Kapal Phinisi 3D");
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f); 
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    // Atur cahaya ambien global
    GLfloat ambient_light[] = {0.3f, 0.3f, 0.3f, 1.0f};
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient_light);
    
    // Atur properti specular untuk cahaya
    GLfloat light_specular[] = {1.0f, 1.0f, 1.0f, 1.0f};
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
    
    // Generate semua geometri
    generateHullMesh();
    generateDeckhouseMesh();
    generateMastMesh({0.0f, 7.5f, 5.0f}, 38.0f, 0.5f, 16); // Tiang utama (lebih depan & lebih tinggi)
    generateMastMesh({0.0f, 10.0f, -8.0f}, 22.0f, 0.4f, 16); // Tiang kedua (lebih belakang & lebih pendek)
    generateRailings();
    generateBowsprit();
    // ---- TAMBAHKAN BLOK INI UNTUK MEMBUAT LAYAR UTAMA ----
    // --- TIANG DEPAN (MAIN MAST) ---
    // --- TIANG DEPAN ---
    glm::vec3 main_p0 = {0.0f, 15.0f, 5.0f};  // Titik bawah di tiang
    glm::vec3 main_p1 = {0.0f, 13.5f, -2.0f}; // Titik bawah di ujung spar
    glm::vec3 main_p2 = {0.0f, 39.0f, -2.0f}; // Titik atas di ujung spar
    glm::vec3 main_p3 = {0.0f, 42.0f, 5.0f};   // Titik ikat atas di tiang
        generateSailMesh(main_p0, main_p1, main_p2, main_p3, 15, sail_vertices, sail_indices);
        generateRailSegment(main_p0, main_p1, 0.15f); // Spar bawah (boom)
        generateRailSegment(main_p3, main_p2, 0.15f); // Spar atas (gaff)
    // Layar Atas Depan
    glm::vec3 main_top_p0 = {0.0f, 30.2f, 5.0f};
    glm::vec3 main_top_p1 = {0.0f, 26.2f, -2.0f};
    glm::vec3 main_top_p2 = {0.0f, 30.0f, 5.0f};
    generateSailMesh(main_top_p0, main_top_p1, main_top_p2, main_top_p2, 10, sail_vertices, sail_indices);

    // --- TIANG BELAKANG (AFT MAST) ---
    // Layar Utama Belakang
    glm::vec3 aft_p0 = {0.0f, 18.0f, -8.0f}; // Y dinaikkan
    glm::vec3 aft_p1 = {0.0f, 15.5f, -15.0f}; // Y dinaikkan
    glm::vec3 aft_p2 = {0.0f, 25.0f, -15.0f};
    glm::vec3 aft_p3 = {0.0f, 30.0f, -8.0f};
    generateSailMesh(aft_p0, aft_p1, aft_p2, aft_p3, 15, sail_vertices, sail_indices);
    generateRailSegment(aft_p0, aft_p1, 0.15f); // Boom bawah
    generateRailSegment(aft_p3, aft_p2, 0.15f); // Gaff atas
    // --- TIGA LAYAR DEPAN (JIBS) - VERSI FINAL YANG TERPISAH ---
    // Tentukan titik-titik utama
    glm::vec3 mast_top   = {0.0f, 35.0f, 5.0f}; // Puncak tiang utama
    glm::vec3 bowsprit_tip = bowsprit_end_point; // Ujung bowsprit

    // Buat titik-titik imajiner di sepanjang tali forestay
    glm::vec3 forestay_p1 = glm::mix(mast_top, bowsprit_tip, 0.4f); // 40% dari atas
    glm::vec3 forestay_p2 = glm::mix(mast_top, bowsprit_tip, 0.7f); // 70% dari atas

    // Tentukan titik-titik ikat di tiang
    glm::vec3 mast_p1 = {0.0f, 28.0f, 5.0f};
    glm::vec3 mast_p2 = {0.0f, 21.0f, 5.0f};
    glm::vec3 mast_p3 = {0.0f, 14.0f, 5.0f};

    // Jib 1 (Paling Atas) - Menempel ke titik forestay_p1
    generateSailMesh(mast_top, mast_p1, forestay_p1, forestay_p1, 1,sail_vertices, sail_indices);

    // Jib 2 (Tengah) - Menempel ke titik forestay_p2
    generateSailMesh(mast_p1, mast_p2, forestay_p2, forestay_p2, 1, sail_vertices, sail_indices);

    // Jib 3 (Bawah) - Menempel ke ujung bowsprit
    generateSailMesh(mast_p2, mast_p3,bowsprit_tip, bowsprit_tip, 1, sail_vertices, sail_indices);
    generateRailSegment(mast_top, bowsprit_tip, 0.08f); // Rail forestay utama
    generateRailSegment(mast_p1, forestay_p1, 0.08f);   // Rail jib 1
    generateRailSegment(mast_p2, forestay_p2, 0.08f);   // Rail jib 2
    generateRailSegment(mast_p3, bowsprit_tip, 0.08f);  // Rail jib 3
    
    generateRigging();
    generateLadders();
    generateCabinLadder();
    

    recalculateAllNormals(); 
    recalculateSailNormals();
    glutDisplayFunc(display);

    exportToObj("F:\\CODING\\CC++\\kapalPhinisi\\kapal_phinisi.obj");
    glutKeyboardFunc(keyboard);
    glutMouseFunc(mouseButton);
    glutMotionFunc(mouseMove);
    glutMouseWheelFunc(mouseWheel);
    glutMainLoop();
    
    return 0;
}