from settings import *
from bird import Bird
from pipe import Pipe
from button import Button
from sys import exit
import random

class Main:
	def __init__(self):
		pygame.init()
		self.screen = pygame.display.set_mode((screen_width, screen_height))
		pygame.display.set_caption(name_game)
		self.clock = pygame.time.Clock()

		self.backgroud = pygame.image.load('img/bg.png')
		self.ground = pygame.image.load('img/ground.png')

		self.bird_group = pygame.sprite.Group()
		self.pipe_group = pygame.sprite.Group()

		self.font = pygame.font.SysFont('Bauhaus 93', 60)

	def draw_text(self, text, font, text_col, x, y):
		img = font.render(text, True, text_col)
		self.screen.blit(img, (x, y))

	def reset_game(self, flappy):
		self.pipe_group.empty()
		flappy.rect.x = 100
		flappy.rect.y = int(screen_height / 2)
		score = 0
		return score

	def run(self):
		ground_scroll = 0
		flappy = Bird(100, int(screen_height / 2))
		self.bird_group.add(flappy)
		button = Button(screen_width//2 - 50, screen_height//2 - 100)
		game_over = False
		last_pipe = pygame.time.get_ticks() - pipe_frequency
		score = 0
		pass_pipe = False
		while True:
			self.clock.tick(fps)
			self.screen.blit(self.backgroud, (0, 0))


			self.bird_group.draw(self.screen)
			self.bird_group.update(game_over)

			self.pipe_group.draw(self.screen)

			self.screen.blit(self.ground, (ground_scroll, ground_cord))

			if len(self.pipe_group) > 0:
				if self.bird_group.sprites()[0].rect.left > self.pipe_group.sprites()[0].rect.left\
					and self.bird_group.sprites()[0].rect.right < self.pipe_group.sprites()[0].rect.right\
					and pass_pipe == False:
					pass_pipe = True
				if pass_pipe == True:
					if self.bird_group.sprites()[0].rect.left > self.pipe_group.sprites()[0].rect.right:
						score += 1
						pass_pipe = False
			
			self.draw_text(str(score), self.font, white, screen_width//2, 20)

			if pygame.sprite.groupcollide(self.bird_group, self.pipe_group, False, False) or flappy.rect.top < 0:
				game_over = True

			if flappy.rect.bottom > 768:
				game_over = True

			if not game_over:

				time_now = pygame.time.get_ticks()
				if time_now - last_pipe > pipe_frequency:
					pipe_height = random.randint(-100, 100)
					btm_pipe = Pipe(screen_width, int(screen_height / 2) + pipe_height, -1)
					top_pipe = Pipe(screen_width, int(screen_height / 2) + pipe_height, 1)
					self.pipe_group.add(btm_pipe)
					self.pipe_group.add(top_pipe)
					last_pipe = time_now



				ground_scroll -= scroll_speed
				if abs(ground_scroll) > 35:
					ground_scroll = 0
				self.pipe_group.update()

			if game_over == True:
				if button.draw(self.screen) == True:
					game_over = False
					score = self.reset_game(flappy)

			for event in pygame.event.get():
				if event.type == pygame.QUIT:
					pygame.quit()
					exit()

			pygame.display.update()


if __name__ == "__main__":
	main = Main()
	main.run()