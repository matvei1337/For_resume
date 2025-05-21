from settings import *

class Pipe(pygame.sprite.Sprite):
	def __init__(self, x, y, positioin):
		pygame.sprite.Sprite.__init__(self)
		self.image = pygame.image.load('img/pipe.png')
		self.rect = self.image.get_rect()

		if positioin == 1:
			self.image = pygame.transform.flip(self.image, False, True)
			self.rect.bottomleft = [x, y - pipe_gap//2]
		else:
			self.rect.topleft = [x, y + pipe_gap//2]

	def update(self):
		self.rect.x -= scroll_speed
		if self.rect.right < 0:
			self.kill()