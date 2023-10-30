static struct sigaction new_sa, old_sa;

new_sa.sa_handler = my_handler;
sigemptyset(&new_handler.sa_mask);

if (sigaction(signo, &new_sa, &old_sa) == -1) {
    /* handle sigaction error */
}