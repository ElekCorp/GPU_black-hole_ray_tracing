#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "black_hole.cuh"
#include "cuda_ray.cuh"

#include "szinsaver.h"

#include <iostream>
#include <fstream>

//#include <windows.h>

#include <chrono>

#include <vector>

#include <cmath>

#include <SDL.h>
#include <SDL_ttf.h>

template <class FP>
int8_t* makeframe(int SZELES, int MAGAS, FP* x, FP* Omega, FP a, FP Q, FP rs, FP errormax, FP de0, FP kepernyo_high, FP kepernyo_tav, FP sugar_ki, FP gyuru_sugar_kicsi, FP gyuru_sugar_nagy, int SZELESregi, int MAGASregi, int ikezd, int jkezd, int iveg);

template <class FP>
FP* makeframe_T(int SZELES, int MAGAS, FP* x, FP* Omega, FP a, FP Q, FP rs, FP errormax, FP de0, FP kepernyo_high, FP kepernyo_tav, FP sugar_ki, FP gyuru_sugar_kicsi, FP gyuru_sugar_nagy, int SZELESregi, int MAGASregi, int ikezd, int jkezd, int iveg);

void console_kep(int MAGAS, int SZELES, int8_t* SZIN);
void device_info(void);
void SDL_szinezo(SDL_Renderer* ren, int8_t* SZIN, int i, int j, int ki, int kj, int MAGAS);
//void fps_kiiro(int fps);

int n_oszto(int SZELES, int MAGAS, int kepernyoSZELES, int kepernyoMAGAS, int n);


int main(int argc, char* argv[])
{

    device_info();

    auto start = std::chrono::high_resolution_clock::now();

    int kepernyoSZELES = 1024;
    int kepernyoMAGAS = 512;
    
    int SZELES = 64;
    int MAGAS = 32;

    if (kepernyoSZELES * MAGAS != kepernyoMAGAS * SZELES)
    {
        std::cout << "renderelendo kep es a kepernyo aranya nem azonos\n";
        return 1;
    }

    float errormax = 0.001f;
    float de0 = 0.01f;

    float rs = 0.05f;//2*rs=m

    float delta_a = 0.0001f;

    float a = 0.0f;//rs / 2 -delta_a;//0.0;
    float Q = 0.0f;



    float t_0 = 0.0f;
    float r_0 = 1.0f;
    float theta_0 = 1.57f + 0.06f;//ne legyen nulla
    float phi_0 = 0.0f;

    float kepernyo_high = 0.5f;
    float kepernyo_tav = 0.75f;//0.4;//0.75

    float sugar_ki = 1.01f;

    float gyuru_sugar_kicsi = 0.1f;
    float gyuru_sugar_nagy = 0.5f;

    float x[D] = { t_0,r_0,theta_0,phi_0 };
    float pi_cucc = float(asin(1) * 2);
    float Omega[D - 1] = { 0,pi_cucc,0 };
    
    


    int SZELESregi = SZELES;
    int MAGASregi = MAGAS;
    
    int ikezd = 0;
    int jkezd = 0;
    int iveg = SZELES;
    

    if (SDL_Init(SDL_INIT_EVERYTHING) < 0) {
        SDL_Log("Nem indithato az SDL: %s", SDL_GetError());
        exit(1);
    }

    SDL_Window* window = SDL_CreateWindow("SDL black hole", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, kepernyoSZELES, kepernyoMAGAS, 0);
    //SDL_Window* window = SDL_CreateWindow("SDL black hole", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, kepernyoSZELES, kepernyoMAGAS, SDL_WINDOW_RESIZABLE);//ezzel at meretezheto az ablak merete de a kirenderelt pixelek nem valtoznak

    if (window == NULL) {
        SDL_Log("Nem hozhato letre az ablak: %s", SDL_GetError());
        exit(1);
    }

    //icon teszt

    SDL_Surface* surface;     // Declare an SDL_Surface to be filled in with pixel data from an image file
    /*Uint16 pixels[16 * 16] = {  //SDL icon
      0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff,
      0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff,
      0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff,
      0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff,
      0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff,
      0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff,
      0x0fff, 0x0aab, 0x0789, 0x0bcc, 0x0eee, 0x09aa, 0x099a, 0x0ddd,
      0x0fff, 0x0eee, 0x0899, 0x0fff, 0x0fff, 0x1fff, 0x0dde, 0x0dee,
      0x0fff, 0xabbc, 0xf779, 0x8cdd, 0x3fff, 0x9bbc, 0xaaab, 0x6fff,
      0x0fff, 0x3fff, 0xbaab, 0x0fff, 0x0fff, 0x6689, 0x6fff, 0x0dee,
      0xe678, 0xf134, 0x8abb, 0xf235, 0xf678, 0xf013, 0xf568, 0xf001,
      0xd889, 0x7abc, 0xf001, 0x0fff, 0x0fff, 0x0bcc, 0x9124, 0x5fff,
      0xf124, 0xf356, 0x3eee, 0x0fff, 0x7bbc, 0xf124, 0x0789, 0x2fff,
      0xf002, 0xd789, 0xf024, 0x0fff, 0x0fff, 0x0002, 0x0134, 0xd79a,
      0x1fff, 0xf023, 0xf000, 0xf124, 0xc99a, 0xf024, 0x0567, 0x0fff,
      0xf002, 0xe678, 0xf013, 0x0fff, 0x0ddd, 0x0fff, 0x0fff, 0xb689,
      0x8abb, 0x0fff, 0x0fff, 0xf001, 0xf235, 0xf013, 0x0fff, 0xd789,
      0xf002, 0x9899, 0xf001, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0xe789,
      0xf023, 0xf000, 0xf001, 0xe456, 0x8bcc, 0xf013, 0xf002, 0xf012,
      0x1767, 0x5aaa, 0xf013, 0xf001, 0xf000, 0x0fff, 0x7fff, 0xf124,
      0x0fff, 0x089a, 0x0578, 0x0fff, 0x089a, 0x0013, 0x0245, 0x0eff,
      0x0223, 0x0dde, 0x0135, 0x0789, 0x0ddd, 0xbbbc, 0xf346, 0x0467,
      0x0fff, 0x4eee, 0x3ddd, 0x0edd, 0x0dee, 0x0fff, 0x0fff, 0x0dee,
      0x0def, 0x08ab, 0x0fff, 0x7fff, 0xfabc, 0xf356, 0x0457, 0x0467,
      0x0fff, 0x0bcd, 0x4bde, 0x9bcc, 0x8dee, 0x8eff, 0x8fff, 0x9fff,
      0xadee, 0xeccd, 0xf689, 0xc357, 0x2356, 0x0356, 0x0467, 0x0467,
      0x0fff, 0x0ccd, 0x0bdd, 0x0cdd, 0x0aaa, 0x2234, 0x4135, 0x4346,
      0x5356, 0x2246, 0x0346, 0x0356, 0x0467, 0x0356, 0x0467, 0x0467,
      0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff,
      0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff,
      0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff,
      0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff, 0x0fff
    };
    surface = SDL_CreateRGBSurfaceFrom(pixels, 16, 16, 16, 16 * 2, 0x0f00, 0x00f0, 0x000f, 0xf000);*/
    char* icon_helye = "icon_black_hole.bmp";
    surface = SDL_LoadBMP(icon_helye);
    // The icon is attached to the window pointer
    SDL_SetWindowIcon(window, surface);

    // ...and the surface containing the icon pixel data is no longer required.
    SDL_FreeSurface(surface);












    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_SOFTWARE);
    if (renderer == NULL) {
        SDL_Log("Nem hozhato letre a megjelenito: %s", SDL_GetError());
        exit(1);
    }
    SDL_RenderClear(renderer);


    /*TTF_Init();
    TTF_Font* font = TTF_OpenFont("AGENCYR.ttf", 24);
    if (!font) {
        SDL_Log("Nem sikerult megnyitni a fontot! %s\n", TTF_GetError());
        exit(1);
    }

    SDL_Surface* felirat;
    SDL_Texture* felirat_t;
    SDL_Rect hova = { 0, 0, 0, 0 };
    SDL_Color szoveg_color = { 255, 255, 255 };*/
    /* ha sajat kodban hasznalod, csinalj belole fuggvenyt! */
    /*felirat = TTF_RenderUTF8_Solid(font, "Teszt123", szoveg_color);
    felirat_t = SDL_CreateTextureFromSurface(renderer, felirat);
    hova.x = 0;
    hova.y = 0;
    hova.w = felirat->w;
    hova.h = felirat->h;
    SDL_RenderCopy(renderer, felirat_t, NULL, &hova);
    */
    
    
    bool quit = false;
    bool lent = false;

    bool kep_frissito = true;

    int mousex;
    int mousey;

    int szamolo = 0;

    std::vector<int> x_coord, y_coord;

    while (!quit)
    {

        Uint64 start_perf = SDL_GetPerformanceCounter();

        SDL_Event event;
        SDL_WaitEvent(&event);//ezzel akkor lép tovább ha kap eventet
        const Uint8* state = SDL_GetKeyboardState(NULL);

        int8_t* SZIN = NULL;

        if ((kepernyoSZELES * MAGAS != kepernyoMAGAS * SZELES) || (kepernyoMAGAS % MAGAS != 0))
        {
            std::cout << "renderelendo kep es a kepernyo aranya nem azonos\n";
            //exit(EXIT_FAILURE);
            SZELES = kepernyoSZELES;
            MAGAS = kepernyoMAGAS;
            
        }
        
        if (kep_frissito == true)
        {
            SZIN = makeframe(SZELES, MAGAS, x, Omega, a, Q, rs, errormax, de0, kepernyo_high, kepernyo_tav, sugar_ki, gyuru_sugar_kicsi, gyuru_sugar_nagy, SZELESregi, MAGASregi, ikezd, jkezd, iveg);
            SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
            SDL_RenderClear(renderer);

            kep_frissito = false;
            //std::cout << "new frame\n";



            for (int ki = 0; ki < SZELES; ++ki)
            {
                for (int kj = 0; kj < MAGAS; ++kj)
                {
                    if (kepernyoSZELES == SZELES && kepernyoMAGAS == MAGAS)
                    {
                        int i = ki;
                        int j = kj;

                        SDL_szinezo(renderer, SZIN, i, j, ki, kj, MAGAS);
                    }
                    else if (kepernyoSZELES % SZELES == 0)
                    {
                        int num = kepernyoSZELES / SZELES;
                        for (int k = 0; k < num; ++k)
                        {
                            int i = num * ki + k;
                            int j = num * kj + k;

                            SDL_szinezo(renderer, SZIN, i, j, ki, kj, MAGAS);
                        }



                    }


                }
            }
            free(SZIN);

            //SDL_RenderCopy(renderer, felirat_t, NULL, &hova);
            SDL_RenderPresent(renderer);
        }

        float addx = 0.01f;

        switch (event.type)
        {
        case SDL_KEYDOWN:
            if (state[SDL_SCANCODE_UP])
            {
                
                x[1] -= addx;
                sugar_ki -= addx;
                
            }
            if (state[SDL_SCANCODE_DOWN])
            {
                x[1] += addx;
                sugar_ki += addx;
                
            }
            if (state[SDL_SCANCODE_W])
            {
                x[2] += addx;
                
            }
            if (state[SDL_SCANCODE_S])
            {
                x[2] -= addx;
                
            }
            if (state[SDL_SCANCODE_LEFT])
            {
                a -= delta_a;
            }
            if (state[SDL_SCANCODE_RIGHT])
            {
                a += delta_a;
                
            }
            if (state[SDL_SCANCODE_X])
            {
                rs -= addx/10;
                
            }
            if (state[SDL_SCANCODE_C])
            {
                rs += addx/10;
                
            }
            if (state[SDL_SCANCODE_D])
            {
                Omega[1] += addx;
                //std::cout << Omega[1] / pi_cucc * 180.0 << "\n";
                
            }
            if (state[SDL_SCANCODE_A])
            {
                Omega[1] -= addx;
                //std::cout << Omega[1] / pi_cucc * 180.0 << "\n";
                
            }
            if (state[SDL_SCANCODE_KP_PLUS])
            {
                errormax *= 1.01f;
                std::cout << errormax << "\n";
                
            }
            if (state[SDL_SCANCODE_KP_MINUS])
            {
                errormax /= 1.01f;
                
            }
            if (state[SDL_SCANCODE_KP_0])
            {
                
                std::cout << "reset frame\n";
                szamolo = 0;

                SZELES = kepernyoSZELES;
                MAGAS = kepernyoMAGAS;

                SZELESregi = kepernyoSZELES;
                MAGASregi = kepernyoMAGAS;


                ikezd = 0;
                jkezd = 0;
                iveg = SZELES;
            }
            if (state[SDL_SCANCODE_KP_1])
            {
                SZELES = kepernyoSZELES;
                MAGAS = kepernyoMAGAS;


                
                //SDL_RenderClear(renderer);
            }
            if (state[SDL_SCANCODE_KP_2])
            {
                

                int oszto = n_oszto(SZELES, MAGAS, kepernyoSZELES, kepernyoMAGAS, 2);
                SZELES = kepernyoSZELES / oszto;
                MAGAS = kepernyoMAGAS / oszto;


                //SDL_RenderClear(renderer);
            }
            if (state[SDL_SCANCODE_KP_3])
            {
                

                int oszto = n_oszto(SZELES, MAGAS, kepernyoSZELES, kepernyoMAGAS, 3);
                SZELES = kepernyoSZELES / oszto;
                MAGAS = kepernyoMAGAS / oszto;
;

                //SDL_RenderClear(renderer);
            }
            if (state[SDL_SCANCODE_KP_4])
            {
                

                int oszto = n_oszto(SZELES, MAGAS, kepernyoSZELES, kepernyoMAGAS, 4);
                SZELES = kepernyoSZELES / oszto;
                MAGAS = kepernyoMAGAS / oszto;

                //SDL_RenderClear(renderer);
            }
            if (state[SDL_SCANCODE_KP_5])
            {
                

                int oszto = n_oszto(SZELES, MAGAS, kepernyoSZELES, kepernyoMAGAS, 5);
                SZELES = kepernyoSZELES / oszto;
                MAGAS = kepernyoMAGAS / oszto;

                //SDL_RenderClear(renderer);
            }
            if (state[SDL_SCANCODE_KP_6])
            {
                

                int oszto = n_oszto(SZELES, MAGAS, kepernyoSZELES, kepernyoMAGAS, 6);
                SZELES = kepernyoSZELES / oszto;
                MAGAS = kepernyoMAGAS / oszto;

                //SDL_RenderClear(renderer);
            }
            if (state[SDL_SCANCODE_KP_7])
            {
                

                int oszto = n_oszto(SZELES, MAGAS, kepernyoSZELES, kepernyoMAGAS, 7);
                SZELES = kepernyoSZELES / oszto;
                MAGAS = kepernyoMAGAS / oszto;

                //SDL_RenderClear(renderer);
            }
            if (state[SDL_SCANCODE_KP_8])
            {
                

                int oszto = n_oszto(SZELES, MAGAS, kepernyoSZELES, kepernyoMAGAS, 8);
                SZELES = kepernyoSZELES / oszto;
                MAGAS = kepernyoMAGAS / oszto;
            
                //SDL_RenderClear(renderer);

            }
            if (state[SDL_SCANCODE_KP_9])
            {
                

                int oszto = n_oszto(SZELES, MAGAS, kepernyoSZELES, kepernyoMAGAS, 9);
                SZELES = kepernyoSZELES / oszto;
                MAGAS = kepernyoMAGAS / oszto;

                //SDL_RenderClear(renderer);

            }
            if (state[SDL_SCANCODE_RETURN])//kep mentese ENTER lenyomassal
            {
                std::cout << "kep_mentese 4k-ban\n";

                SZELES = kepernyoSZELES*2;
                MAGAS = kepernyoMAGAS*2;
                float error_szabalyzo = 10;
                errormax /= error_szabalyzo;

                std::cout << "press D to double\nor any for float\n";
                SDL_Event event_double;
                do
                {
                    SDL_WaitEvent(&event_double);
                } while (state[SDL_SCANCODE_RETURN]);
                
                SDL_WaitEvent(&event_double);

                switch (event_double.type)
                {
                case SDL_KEYDOWN:
                    if (state[SDL_SCANCODE_D])
                    {
                        std::cout << "double kep\n";
                        double errormax_d = double(errormax);
                        double de0_d = double(de0);

                        double rs_d = double(rs);//2*rs=m

                        double delta_a_d = double(delta_a);

                        double a_d = double(a);//rs / 2 -delta_a;//0.0;
                        double Q_d = double(Q);

                        double kepernyo_high_d = double(kepernyo_high);
                        double kepernyo_tav_d = double(kepernyo_tav);//0.4;//0.75

                        double sugar_ki_d = double(sugar_ki);

                        double gyuru_sugar_kicsi_d = double(gyuru_sugar_kicsi);
                        double gyuru_sugar_nagy_d = double(gyuru_sugar_nagy);

                        double x_d[D] = { double(x[0]),double(x[1]),double(x[2]),double(x[3]) };
                        double Omega_d[D - 1] = { double(Omega[0]),double(Omega[1]),double(Omega[2]) };

                        double* SZIN_fp = makeframe_T(SZELES, MAGAS, x_d, Omega_d, a_d, Q_d, rs_d, errormax_d, de0_d, kepernyo_high_d, kepernyo_tav_d, sugar_ki_d, gyuru_sugar_kicsi_d, gyuru_sugar_nagy_d, SZELESregi, MAGASregi, ikezd, jkezd, iveg);

                        datasaver_T(SZIN_fp, SZELES, MAGAS, "kep_double.dat");
                        free(SZIN_fp);
                        std::cout << "kesz a kep_double.dat\n";
                    }
                    else
                    {
                        float* SZIN_fp = makeframe_T(SZELES, MAGAS, x, Omega, a, Q, rs, errormax, de0, kepernyo_high, kepernyo_tav, sugar_ki, gyuru_sugar_kicsi, gyuru_sugar_nagy, SZELESregi, MAGASregi, ikezd, jkezd, iveg);

                        datasaver_T(SZIN_fp, SZELES, MAGAS, "kep.dat");
                        free(SZIN_fp);
                        std::cout << "kesz a kep.dat\n";
                    }
                    
                    break;

                default:
                    float* SZIN_fp = makeframe_T(SZELES, MAGAS, x, Omega, a, Q, rs, errormax, de0, kepernyo_high, kepernyo_tav, sugar_ki, gyuru_sugar_kicsi, gyuru_sugar_nagy, SZELESregi, MAGASregi, ikezd, jkezd, iveg);

                    datasaver_T(SZIN_fp, SZELES, MAGAS, "kep.dat");
                    free(SZIN_fp);
                    std::cout << "kesz a kep.dat\n";
                    break;
                }


                

                
                errormax *= error_szabalyzo;
                
                int oszto = n_oszto(SZELES, MAGAS, kepernyoSZELES, kepernyoMAGAS, 3);
                SZELES = kepernyoSZELES / oszto;
                MAGAS = kepernyoMAGAS / oszto;



            }
            if (state[SDL_SCANCODE_ESCAPE])
            {
                quit = true;
            }

            kep_frissito = true;

            break;

        case SDL_MOUSEBUTTONDOWN:
            if (lent == false)
            {
                lent = true;
            }
            
            break;

        case SDL_MOUSEBUTTONUP:
            

            if (lent == true)
            {
                SDL_GetMouseState(&mousex, &mousey);
                x_coord.push_back(mousex);
                y_coord.push_back(mousey);
                std::cout << mousex << "\t" << mousey << "\n";


                ++szamolo;
                if (szamolo == 2)
                {
                    size_t indmax = x_coord.size();

                    SZELESregi = kepernyoSZELES;
                    MAGASregi = kepernyoMAGAS;


                    ikezd = min(x_coord[indmax - 1], x_coord[indmax - 2]);
                    jkezd = min(y_coord[indmax - 1], y_coord[indmax - 2]);
                    iveg = max(x_coord[indmax - 1], x_coord[indmax - 2]);

                    kep_frissito = true;
                }
                else if (szamolo % 2 == 0)
                {
                    size_t indmax = x_coord.size();

                    double alpha = double(kepernyoSZELES) / (double(iveg) - double(ikezd));

                    //std::cout << "alpha:" << alpha << "\n";

                    SZELESregi *= alpha;
                    MAGASregi *= alpha;

                    ikezd *= alpha;
                    jkezd *= alpha;
                    iveg *= alpha;
                    //std::cout << "cucc:" << double(kepernyoSZELES) + (double(ikezd) - double(iveg)) << "\n";
                    std::cout << SZELESregi << "\n";

                    ikezd += min(x_coord[indmax - 1], x_coord[indmax - 2]);
                    jkezd += min(y_coord[indmax - 1], y_coord[indmax - 2]);
                    iveg = ikezd + max(x_coord[indmax - 1], x_coord[indmax - 2]);

                    kep_frissito = true;
                }

                


                lent = false;
            }

            
            
            

            break;
            

            break;
        case SDL_QUIT:
            quit = true;
            break;
        }


        Uint64 end_perf = SDL_GetPerformanceCounter();

        float elapsedMS = (end_perf - start_perf) / (float)SDL_GetPerformanceFrequency();// *1000.0f;
        int fps = int(1.0f / elapsedMS);
        const char fps_szoveg[] = { char(48 + fps / 100),char(48 + (fps % 100) / 10),char(48 + fps % 10) };


        /*SDL_FreeSurface(felirat);
        SDL_DestroyTexture(felirat_t);

        felirat = TTF_RenderUTF8_Solid(font, fps_szoveg, szoveg_color);
        felirat_t = SDL_CreateTextureFromSurface(renderer, felirat);
        */
        //SDL_RenderCopy(renderer, felirat_t, NULL, &hova);
        //SDL_RenderPresent(renderer);

        //std::cout << "Current FPS: " << (1.0f / elapsedMS) << std::endl;
    }

    /*SDL_FreeSurface(felirat);
    SDL_DestroyTexture(felirat_t);*/

    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}

template <class FP>
int8_t* makeframe(int SZELES, int MAGAS, FP* x, FP* Omega, FP a, FP Q, FP rs, FP errormax, FP de0, FP kepernyo_high, FP kepernyo_tav, FP sugar_ki, FP gyuru_sugar_kicsi, FP gyuru_sugar_nagy, int SZELESregi, int MAGASregi, int ikezd, int jkezd, int iveg)
{

    FP* x_d = NULL;
    FP* Omega_d = NULL;
    cudaError_t err = cudaSuccess;

    //auto start = std::chrono::high_resolution_clock::now();

    err = cudaMalloc(&x_d, D * sizeof(FP));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "\n cudaMalloc hiba x_d foglalasakor %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaMalloc(&Omega_d, (D - 1) * sizeof(FP));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "\n cudaMalloc hiba Omega_d foglalasakor %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(x_d, x, D * sizeof(FP), cudaMemcpyHostToDevice);
    err = cudaMemcpy(Omega_d, Omega, (D - 1) * sizeof(FP), cudaMemcpyHostToDevice);

    int8_t* SZIN_d = NULL;

    err = cudaMalloc(&SZIN_d, size_t(SZELES) * size_t(MAGAS) * sizeof(int8_t));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "\n cudaMalloc hiba SZIN_d foglalasakor %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    dim3 threadsPerBlock(16, 16);
    int xdim = SZELES / threadsPerBlock.x, ydim = MAGAS / threadsPerBlock.y;

    if (SZELES % threadsPerBlock.x != 0)
    {
        xdim = SZELES / threadsPerBlock.x + 1;
    }
    if (MAGAS % threadsPerBlock.y != 0)
    {
        ydim = MAGAS / threadsPerBlock.y + 1;
    }

    dim3 numBlocks(xdim, ydim);

    ray_step <<<numBlocks, threadsPerBlock >>> (SZIN_d, SZELES, MAGAS, x_d, Omega_d, a, Q, rs, errormax, de0, kepernyo_high, kepernyo_tav, sugar_ki, gyuru_sugar_kicsi, gyuru_sugar_nagy, SZELESregi, MAGASregi, ikezd, jkezd, iveg);

    cudaDeviceSynchronize();

    err = cudaPeekAtLastError();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "\nkernel:%s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    int8_t* SZIN = NULL;
    SZIN = (int8_t*)malloc(SZELES * MAGAS * sizeof(int8_t));

    err = cudaMemcpy(SZIN, SZIN_d, SZELES * MAGAS * sizeof(int8_t), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "\ncudaMemcpy SZIN,SZIN_d, DeviceToHost nal:%s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    cudaFree(SZIN_d);
    cudaFree(x_d);
    cudaFree(Omega_d);

    //auto end = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    //std::cout << "Teljes lefutasi ido\n";
    //std::cout << "\n\n\n\n\n\n\n";
    //std::cout << "\n\n\n\n\n\n\n";
    //std::cout << "\n\n\n\n\n\n\n";
    //std::cout << "\n\n\n\n\n\n\n";
    
    //std::cout << double(duration.count()) / 1000000 << "sec\n"<<(1/(double(duration.count()) / 1000000))<<"fps\n";

    //console_kep(MAGAS, SZELES, SZIN);

    return SZIN;

}


template <class FP>
FP* makeframe_T(int SZELES, int MAGAS, FP* x, FP* Omega, FP a, FP Q, FP rs, FP errormax, FP de0, FP kepernyo_high, FP kepernyo_tav, FP sugar_ki, FP gyuru_sugar_kicsi, FP gyuru_sugar_nagy, int SZELESregi, int MAGASregi, int ikezd, int jkezd, int iveg)//ekkor a SZIN egy FP* es a homersekletet reprezentalja
{

    FP* x_d = NULL;
    FP* Omega_d = NULL;
    cudaError_t err = cudaSuccess;

    //auto start = std::chrono::high_resolution_clock::now();

    err = cudaMalloc(&x_d, D * sizeof(FP));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "\n cudaMalloc hiba x_d foglalasakor %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaMalloc(&Omega_d, (D - 1) * sizeof(FP));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "\n cudaMalloc hiba Omega_d foglalasakor %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    err = cudaMemcpy(x_d, x, D * sizeof(FP), cudaMemcpyHostToDevice);
    err = cudaMemcpy(Omega_d, Omega, (D - 1) * sizeof(FP), cudaMemcpyHostToDevice);

   FP* SZIN_d = NULL;

    err = cudaMalloc(&SZIN_d, size_t(SZELES) * size_t(MAGAS) * sizeof(FP));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "\n cudaMalloc hiba SZIN_d foglalasakor %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    dim3 threadsPerBlock(16, 16);
    int xdim = SZELES / threadsPerBlock.x, ydim = MAGAS / threadsPerBlock.y;

    if (SZELES % threadsPerBlock.x != 0)
    {
        xdim = SZELES / threadsPerBlock.x + 1;
    }
    if (MAGAS % threadsPerBlock.y != 0)
    {
        ydim = MAGAS / threadsPerBlock.y + 1;
    }

    dim3 numBlocks(xdim, ydim);

    ray_step_T << <numBlocks, threadsPerBlock >> > (SZIN_d, SZELES, MAGAS, x_d, Omega_d, a, Q, rs, errormax, de0, kepernyo_high, kepernyo_tav, sugar_ki, gyuru_sugar_kicsi, gyuru_sugar_nagy, SZELESregi, MAGASregi, ikezd, jkezd, iveg);

    cudaDeviceSynchronize();

    err = cudaPeekAtLastError();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "\nkernel:%s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    FP* SZIN = NULL;
    SZIN = (FP*)malloc(SZELES * MAGAS * sizeof(FP));

    err = cudaMemcpy(SZIN, SZIN_d, SZELES * MAGAS * sizeof(FP), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "\ncudaMemcpy SZIN,SZIN_d, DeviceToHost nal:%s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    cudaFree(SZIN_d);
    cudaFree(x_d);
    cudaFree(Omega_d);

    //auto end = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    //std::cout << "Teljes lefutasi ido\n";
    //std::cout << "\n\n\n\n\n\n\n";
    //std::cout << "\n\n\n\n\n\n\n";
    //std::cout << "\n\n\n\n\n\n\n";
    //std::cout << "\n\n\n\n\n\n\n";

    //std::cout << double(duration.count()) / 1000000 << "sec\n"<<(1/(double(duration.count()) / 1000000))<<"fps\n";

    //console_kep(MAGAS, SZELES, SZIN);

    return SZIN;

}

/*
void console_kep(int MAGAS,int SZELES, int8_t* SZIN)
{

    HDC hdc = GetDC(GetConsoleWindow());

    COLORREF black = RGB(0, 0, 0);
    COLORREF white = RGB(255, 255, 255);
    COLORREF green = RGB(0, 255, 0);
    COLORREF red = RGB(255, 0, 0);
    COLORREF blue = RGB(0, 0, 255);

    if (SZELES <= 1024 && MAGAS <= 512)
    {
        for (int i = 0; i < SZELES; ++i)
        {
            for (int j = 0; j < MAGAS; ++j)
            {
                if (SZIN[i * MAGAS + j] == 1)
                {
                    SetPixel(hdc, i, j, black);
                }
                else if (SZIN[i * MAGAS + j] == -1)
                {
                    SetPixel(hdc, i, j, white);
                }
                else if (SZIN[i * MAGAS + j] == 0)
                {
                    SetPixel(hdc, i, j, green);
                }
                else if (SZIN[i * MAGAS + j] == 2)
                {
                    SetPixel(hdc, i, j, red);
                }
                else if (SZIN[i * MAGAS + j] == 3)
                {
                    SetPixel(hdc, i, j, blue);
                }
                else
                {

                }
            }
        }
    }
}
*/
void device_info(void)
{
    int num_gpus = 0;
    cudaGetDeviceCount(&num_gpus);

    if (num_gpus == 0)
    {
        std::cout << "no capable GPU\n";
        exit(EXIT_FAILURE);
    }

    //printf("number of host CPUs:\t%d\n", omp_get_num_procs());
    printf("number of CUDA devices:\t%d\n", num_gpus);

    for (int i = 0; i < num_gpus; i++)
    {
        cudaDeviceProp dprop;
        cudaGetDeviceProperties(&dprop, i);
        printf("%d: %s\n", i, dprop.name);
        printf("%d: L2 Cache Size:%d bytes\n", i, dprop.l2CacheSize);


        int driverVersion = 0, runtimeVersion = 0;

        cudaDriverGetVersion(&driverVersion);
        cudaRuntimeGetVersion(&runtimeVersion);
        printf("  CUDA Driver Version / Runtime Version          %d.%d / %d.%d\n",
            driverVersion / 1000, (driverVersion % 100) / 10,
            runtimeVersion / 1000, (runtimeVersion % 100) / 10);
        printf("  CUDA Capability Major/Minor version number:    %d.%d\n",
            dprop.major, dprop.minor);

    }
}

void SDL_szinezo(SDL_Renderer* ren, int8_t* SZIN, int i, int j, int ki, int kj,int MAGAS)
{
    uint8_t temp[3] = { 0, 0, 0 };

    if (SZIN[ki * MAGAS + kj] == 1)
    {
        temp[0] = 0;
        temp[1] = 0;
        temp[2] = 0;
    }
    else if (SZIN[ki * MAGAS + kj] == -1)
    {
        temp[0] = 255;
        temp[1] = 255;
        temp[2] = 255;
    }
    else if (SZIN[ki * MAGAS + kj] == 0)
    {
        temp[0] = 0;
        temp[1] = 255;
        temp[2] = 0;
    }
    else if (SZIN[ki * MAGAS + kj] == 2)
    {
        temp[0] = 255;
        temp[1] = 0;
        temp[2] = 0;
    }
    else if (SZIN[ki * MAGAS + kj] == 3)
    {
        temp[0] = 0;
        temp[1] = 0;
        temp[2] = 255;
    }

    SDL_SetRenderDrawColor(ren, temp[0], temp[1], temp[2], 255);
    SDL_RenderDrawPoint(ren, i, j);
}

int n_oszto(int SZELES, int MAGAS, int kepernyoSZELES, int kepernyoMAGAS, int n)
{
    int oszto = 1;
    int kis_kep=kepernyoMAGAS;
    int num = 1;
    if (kepernyoMAGAS > kepernyoSZELES)
    {
        kis_kep = kepernyoSZELES;
    }
    while (num < n)
    {
        ++oszto;
        if (kis_kep % oszto == 0)
        {
            ++num;
        }

    }

    if ((kepernyoMAGAS % oszto != 0) || (kepernyoSZELES % oszto != 0))
    {
        std::cout << oszto << "A kepernyoMAGAS vagy kepernyoSZELES nem oszthato a viszatert oszto-val\n";
    }
    std::cout << oszto << "\n";
        

    return oszto;
}
