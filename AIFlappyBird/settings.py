import pygame
import random
pygame.font.init()

WIDTH = 600
HEIGHT = 800
FPS = 60

BIRD_IMGS = [pygame.transform.scale2x(pygame.image.load('img/bird1.png')),
	pygame.transform.scale2x(pygame.image.load('img/bird2.png')),
	pygame.transform.scale2x(pygame.image.load('img/bird3.png'))
]
PIPE_IMG = pygame.transform.scale2x(pygame.image.load('img/pipe.png'))
GROUND_IMG = pygame.transform.scale2x(pygame.image.load('img/ground.png'))
BG_IMG = pygame.transform.scale2x(pygame.image.load('img/bg.png'))

STAT_FONT = pygame.font.SysFont('comicsans', 50)