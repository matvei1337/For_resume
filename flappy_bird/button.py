from settings import *

class Button():
	def __init__(self, x, y):
		self.image = pygame.image.load('img/restart.png')
		self.rect = self.image.get_rect()
		self.rect.topleft = (x, y)

	def draw(self, screen):
		action = False

		pos = pygame.mouse.get_pos()

		if self.rect.collidepoint(pos):
			if pygame.mouse.get_pressed()[0] == 1:
				action = True

		screen.blit(self.image, (self.rect.x, self.rect.y))

		return action