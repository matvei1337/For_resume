from settings import *

class Bird(pygame.sprite.Sprite):
	def __init__(self, x, y):
		pygame.sprite.Sprite.__init__(self)
		self.images = []
		self.index = 0
		self.counter = 0
		for num in range(1, 4):
			img = pygame.image.load(f'img/bird{num}.png')
			self.images.append(img)
		self.image = self.images[self.index]
		self.rect = self.image.get_rect()
		self.rect.center = [x, y]
		self.velocity = 0
		self.clicked = False

	def update(self, game_over):
		self.velocity += 0.5
		if self.velocity > 8:
			self.velocity = 8
		if self.rect.bottom < 768:
			self.rect.y += int(self.velocity)


		if pygame.key.get_pressed()[K_SPACE] == 1 and not self.clicked:
			self.clicked = True
			self.velocity = -10

		if pygame.key.get_pressed()[K_SPACE] == 0:
			self.clicked = False

		self.counter += 1
		flap_cooldown = 5

		if self.counter > flap_cooldown:
			self.counter = 0
			self.index = (self.index + 1) % 3

		self.image = self.images[self.index]
		self.image = pygame.transform.rotate(self.images[self.index], -2*self.velocity if not game_over else -90)
